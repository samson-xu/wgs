#!/usr/bin/perl
# See http://factera.stanford.edu/files/license.txt for Usage Agreement and Licensing Details

$|++;

use List::Util qw[min max];
use Statistics::Descriptive;
use File::Spec;
use Symbol qw(gensym);
use Getopt::Std;
use File::Basename;
use IPC::Open3;

#DESCRIPTION===================================================================================================
#vars-fusions-factera.pl Fusion And Chromosomal Translocation Enumeration and Recovery Algorithm
#A tool for rapid identification of genomic fusions in paired-end targeted sequencing data.
#==============================================================================================================

#METHOD========================================================================================================
# 1) Identify improperly paired reads and classify into unique gene/exon pairs
# 2) Find partially mapped reads (soft-clipped) within a reasonable interval of (1)
# 3) Use reads from (2) to determine breakpoint
# 4) BLAST all reads against putative fusion sequences to calculate depth statistics
#==============================================================================================================

#DEPENDENCIES==================================================================================================
# Must have the following executables in path:
# SAMtools, twoBitToFa <http://hgdownload.cse.ucsc.edu/admin/exe/>, blastn, makeblastdb
# Perl dependencies:
# Statistics::Descriptive; Perl 5 core modules = List::Util, IPC::Open3, File::Spec, Symbol, Getopt::Std
#==============================================================================================================

my $version = "1.4.4";

my $size = @ARGV;

my %opts = (r=>5, m=>2, x=>5, s=>1, f=>0.9, S=>0.95, p=>4, k=>12, b=>500, a=>50, c=>16);
my %opts1;
getopts('o:r:m:x:s:f:S:p:veCFk:b:p:a:tc:', \%opts1);

die("
FACTERA version $version by Aaron M. Newman (amnewman\@stanford.edu)

Purpose:
    A tool for detection of fusions in paired-end targeted (or genome-wide) sequencing data.

Usage:
    factera.pl [options] <tumor.bam> <exons.bed> <genome.2bit> <optional: targets.bed>

    <tumor.bam>     Paired reads mapped with BWA or another soft-clip capable aligner.
                    NOTE: Must be indexed for factera to function properly.
    <exons.bed>     Genomic coordinates with gene/exon names in fourth column.
                    Data sources include RefSeq, UCSC, and/or GENCODE.
    <genome.2bit>   Two bit reference sequence (e.g., hg19.2bit).
                    NOTE: twoBitToFa application must be in PATH.
    <targets.bed>   Restrict search to specified genome coordinates.

Options (defaults in parentheses):
    -o <dir>  output directory (tumor.bam directory)
    -r <int>  minimum number of breakpoint-spanning reads required for output ($opts{r})
    -m <int>  minimum number of discordant reads required for each candidate fusion ($opts{m})
    -x <int>  maximum number of breakpoints to examine for each pair of genomic regions ($opts{x})
    -s <int>  minimum number of reads with the same breakpoint ($opts{s})
    -f <0-1>  minimum fraction of read bases required for alignment to fusion template ($opts{f})
    -S <0-1>  minimum similarity required for read to match fusion template ($opts{S})
    -k <int>  k-mer size for fragment comparison ($opts{k} bases)
    -c <int>  minimum size of soft-clipped region to consider ($opts{c} bases)
    -b <int>  number of bases flanking breakpoint for fusion template ($opts{b})
    -p <int>  number of threads for blastn search ($opts{p}; 10 or more recommended)
    -a <int>  number of bases flanking breakpoint to provide in output ($opts{a})
    -e        disable grouping of input coordinates by column 4 of exons.bed (off)
    -v        disable verbose output (off)
    -t        disable running time output (off)
    -C        disable addition of 'chr' prefix to chromosome names (off)
    -F        force remake of BLAST database for a particular input (off)

For detailed help, see http://factera.stanford.edu

") if($size == 0);

my $options='';
foreach my $opt(keys %opts1){
    $options.='|'.$opt.$opts1{$opt};
    $opts{$opt}=$opts1{$opt};
    if($opt eq "o") {
        die("\'$opts1{$opt}\' is not a directory. Abort!\n") if(!(-d $opts1{$opt}));
        print "Output directory: $opts1{$opt}\n";
    }
    if($opt eq "r") {
        die("Minimum breakpoint-spanning reads needed must be an integer. Abort!\n") if is_integer($opts1{$opt}) == 0;
        die("Minimum breakpoint-spanning reads needed must be >0. Abort!\n") if $opts1{$opt}<1;
        print "Minimum breakpoint-spanning reads needed for output set to: $opts1{$opt}\n";
    }
    if($opt eq "m") {
        die("Minimum discordant pairs for candidate fusion must be an integer. Abort!\n") if is_integer($opts1{$opt}) == 0;
        die("Minimum discordant pairs for candidate fusion must be >0. Abort!\n") if $opts1{$opt}<1;
        print "Minimum discordant pairs for candidate fusion set to: $opts1{$opt}\n";
    }
    if($opt eq "x") {
        die("Minimum breakpoints for each gene pair must be an integer. Abort!\n") if is_integer($opts1{$opt}) == 0;
        die("Maximum breakpoints for each gene pair must be >0. Abort!\n") if $opts1{$opt}<1;
        print "Maximum breakpoints for each gene pair set to: $opts1{$opt}\n";
    }
    if($opt eq "s") {
        die("Minimum number of reads with same breakpoint must be an integer. Abort!\n") if is_integer($opts1{$opt}) == 0;
        die("Minimum number of reads with same breakpoint must be >0. Abort!\n") if $opts1{$opt}<1;
        print "Minimum number of reads with same breakpoint set to: $opts1{$opt}\n";
    }
    if($opt eq "f") {
        die("Minimum fraction of matching bases for fusion template must be a fraction (0-1). Abort!\n") if ($opts1{$opt}<0 || $opts1{$opt}>1);
        print "Minimum fraction of matching bases for fusion template set to: $opts1{$opt}\n";
    }
    if($opt eq "S") {
        die("Minimum similarity between reads and fusion template must be a fraction (0-1). Abort!\n") if ($opts1{$opt}<0 || $opts1{$opt}>1);
        print "Minimum similarity between reads and fusion template set to: $opts1{$opt}\n";
    }
    if($opt eq "k") {
        die("k-mer size for fragment comparison must be an integer. Abort!\n") if is_integer($opts1{$opt}) == 0;
        die("k-mer size for fragment comparison must be >0. Abort!\n") if $opts1{$opt}<1;
        print "k-mer size for fragment comparison set to: $opts1{$opt}\n";
    }
    if($opt eq "c") {
        die("Minimum size of soft-clipped read must be an integer. Abort!\n") if is_integer($opts1{$opt}) == 0;
        die("Minimum size of soft-clipped read must be >0. Abort!\n") if $opts1{$opt}<1;
        print "Minimum size of soft-clipped read set to: $opts1{$opt}\n";
    }
    if($opt eq "b") {
        die("Number of bases flanking breakpoint of fusion template must be an integer. Abort!\n") if is_integer($opts1{$opt}) == 0;
        die("Number of bases flanking breakpoint of fusion template must be >0. Abort!\n") if $opts1{$opt}<1;
        print "Number of bases flanking breakpoint of fusion template set to: $opts1{$opt}\n";
    }
    if($opt eq "p") {
        die("Number of threads for blastn search must be an integer. Abort!\n") if is_integer($opts1{$opt}) == 0;
        die("Number of threads for blastn search must be >0. Abort!\n") if $opts1{$opt}<1;
        print "Number of threads for blastn search set to: $opts1{$opt}\n";
    }
    if($opt eq "a") {
        die("Number of bases flanking breakpoint in output must be an integer. Abort!\n") if is_integer($opts1{$opt}) == 0;
        die("Number of bases flanking breakpoint in output must be >0. Abort!\n") if $opts1{$opt}<1;
        print "Number of bases flanking breakpoint in output set to: $opts1{$opt}\n";
    }
    if($opt eq "e") {print "Disable grouping of input coordinates (e.g., allow fusions between exons)\n";}
    if($opt eq "v") {print "Suppress verbose output\n";}
    if($opt eq "t") {print "Suppress running time output\n";}
    if($opt eq "C") {print "Disable addition of \'chr\' prefix to chromosome names\n";}
    if($opt eq "F") {print "Force remake of blast database\n";}
}

#INPUT=========================================================================================================
my $bam; #bam file
my $exon; #exon coordinate bed file [chr <tab> start <tab> end <tab> genename]
my $twobit; #2bit genome file (e.g., hg19.2bit)
my $targets = 0; #restrict analysis to targeted regions (.bed) [use 0 if non-targeted]
#==============================================================================================================

die("Missing input files; Abort!\n") unless (@ARGV > 2);

#get input files
for (my $itor = 0; $itor < @ARGV; $itor++)
{
    my $arg = $ARGV[$itor];
    if($itor == 0) {$bam = $arg;}
    if($itor == 1) {$exon = $arg;}
    if($itor == 2) {$twobit = $arg;}
    if($itor == 3) {
        $targets = $arg;
        print "Restrict results to regions specified in: $targets\n";
    }
}

#PARAMETERS============================================================================================
my $MINSPANNINGREADS = $opts{r}; #minimum number of breakpoint-spanning reads needed to output a fusion
my $MINIMPROPERREADS = $opts{m}; #minimum number of discordant reads needed to consider putative fusion
my $MAXBPS2EXAMINE = $opts{x}; #maximum number of putative breakpoints to consider for each unique gene pair
my $MINBPSUPPORT = $opts{s}; #minimum number of reads spanning putative breakpoint to proceed
my $MINBLASTFRAC = $opts{f}; #fraction of read length to determine minimum alignment length for matching sequences
my $VERBOSE = 1; #if 1, verbose output, including running time; o.w. only results will be printed
if(exists($opts{v})) {$VERBOSE = 0;}
my $TIME = 1; #if 1, print running time
if(exists($opts{t})) {$TIME = 0;}
my $USE_ALL_COORS = 0; #if 0, cluster input coordinates (second argument) by 4th column of bed file (e.g., gene symbol); o.w. treat all coordinates separately (i.e., each line in bed file will be distinct; e.g., to find fusions between exons, rather than between genes)
if(exists($opts{e})) {$USE_ALL_COORS = 1;}
my $makeblastdb = 1; #will set to 0 if a $bam.nhr file is detected
my $FORCEREMAKE = 0; #if 1, force remake of blast database if one already exists
if(exists($opts{F})) {$FORCEREMAKE = 1;}
my $disablechrprefix = 0; #if 1, do not add "chr" prefix to chromosome names
if(exists($opts{C})) {$disablechrprefix = 1;}
my $k = $opts{k}; #k-mer size
my $buffer = $opts{b}; #size of genomic padding around breakpoint
my $pad = $opts{a}; #number of bases flanking breakpoint to provide in output
my $clipsize = $opts{c}; #only keep clipped regions at least 16 bp long by default (only 1 in 4.3B by random chance)
my $minsimilarity = $opts{S}; #minimum %similarity for read to match fusion template
my $BLASTTHREADS = $opts{p}; #number of threads to assign to blastn
my $OUTPUTDIR = dirname($bam); #directory to write output files
my $READLEN = 101; #default read length (will be modified if needed base on input)
my $READLENBIT = 0; #if zero, check input read length and if different, modify READLEN
if(exists($opts{o})) {$OUTPUTDIR = $opts{o};}
die("k-mer size ($k) cannot be greater than minimum size of soft-clipped read ($clipsize). Abort!\n") if $k > $clipsize;
#==============================================================================================================
#write parameters used and arguments to disk
my $params = basename($bam);
$params =~ s/bam/factera.parameters.txt/g;
$params = "$OUTPUTDIR/$params";
open("output0", ">$params") or die $!;
my $v_ = $VERBOSE;
if($v_ == 1) {$v_ = 0;}
my $t_ = $TIME;
if($t_ == 1) {$t_ = 0;}
my $Ta_ = $targets;
if($Ta_ eq "0") {$Ta_ = "null";}
print output0 "Parameter (symbolic argument)\tValue
Input BAM:\t$bam
Genomic coordinates:\t$exon
Reference genome:\t$twobit
Targeted regions:\t$Ta_
Output directory (o):\t$OUTPUTDIR
Minimum breakpoint-spanning reads for each fusion (r):\t$MINSPANNINGREADS
Minimum discordant reads for each candidate fusion (m):\t$MINIMPROPERREADS
Maximum breakpoints to examine for each read pair (x):\t$MAXBPS2EXAMINE
Minimum reads with same breakpoint required (s):\t$MINBPSUPPORT
Minimum fraction of read bases for fusion template alignment (f):\t$MINBLASTFRAC
Minimum % similarity for read to match fusion template (S):\t$minsimilarity
K-mer size for fragment comparison (k):\t$k
Minimum size of soft-clipped read region (c):\t$clipsize
Number of bases flanking breakpoint for fusion template (b):\t$buffer
Number of threads for blastn search (p):\t$BLASTTHREADS
Number of bases flanking breakpoint to provide in output (a):\t$pad
Disable grouping of input coordinates by column 4 of input bed (e):\t$USE_ALL_COORS
Disable verbose output (v):\t$v_
Disable running time output (t):\t$t_
Force remake of blast database (F):\t$FORCEREMAKE
";

close(output0);

#if blast database already exists, don't make it again...

#write reads in fasta format to create blast database of improperly paired, soft-clipped, and unmapped reads
my $blastdb = basename($bam);
$blastdb =~ s/bam/factera.blastreads.fa/g;
$blastdb = "$OUTPUTDIR/$blastdb";
if(-e $blastdb && $FORCEREMAKE == 0){
    $makeblastdb = 0;
    print "Old blast database detected. Use -F argument to force remake.\n";
}

#==============================================================================================================


my @exons_start = (); #store exon start coordinates
my %coors2gene = (); #store gene symbol for every unique exon [format: chr_start_end]
my %chrom = (); #store index for each chromosome
my %getothercoor = (); #get end coordinate if input is start and vice versa [format: gene_start or gene_end]

my $blastitor = 0; #iterate so every subject has unique identifier

my $now = time; #running time
 
open(NULL, ">", File::Spec->devnull); #suppress standard error for mpileup

#read in exon coors: store coordinates in sorted arrays for each chromosome
open(FILE, $exon) or die $!;
while(<FILE>){
    my $line = $_;
    chomp($line);
    my @tokens = split("\t", $line);
    my $chr = addprefix($tokens[0]);
    my $start = $tokens[1];
    my $end = $tokens[2];
    my $gene = $tokens[3];
    if(!exists($chrom{$chr})){
	$size = (keys %chrom);
	$chrom{$chr} = $size;
    }
    my $chrindex = $chrom{$chr};  #retrieve chromosome index
    push(@{$exons_start[$chrindex]}, $start); #add start coordinates
    if(!exists($coors2gene{"$chr $start $end"})) {
        my $insert = $gene;
        if($USE_ALL_COORS == 1) {$insert = "$gene:$start";}
        $coors2gene{"$chr $start $end"} = $insert;
    }
    $getothercoor{"$chr $start"} = $end;
    $getothercoor{"$chr $end"} = $start;
}
if($VERBOSE == 1) {print "Done loading genomic coordinates\n";}

#sort coordinates for each chromosome
for(my $itor = 0; $itor < @exons_start; $itor++){
    @{$exons_start[$itor]} = sort { $a <=> $b } @{$exons_start[$itor]};
}
if($VERBOSE == 1) {print "Done sorting genomic coordinates\n";}

my $impropname = basename($bam);
$impropname =~ s/bam/factera.discordantpair.details.txt/g;
$impropname = "$OUTPUTDIR/$impropname";

open("output", ">$impropname") or die $!;

my $depthname = basename($bam);
$depthname =~ s/bam/factera.discordantpair.depth.txt/g;
$depthname = "$OUTPUTDIR/$depthname";

open(output2, ">$depthname") or die $!;

print output "GENE1\tGENE2\tFUSION\tCHR1\tCHR2\tGENE1_START\tGENE2_START\tDIST\tREAD1_POS\tREAD2_POS\n";

print output2 "FUSION\tDEPTH\n";



if($makeblastdb == 1) {open(blast, ">$blastdb") or die $!;}

#collect all ids of mapped improper pairs in targeted region 
my $insert = "-L $targets";
if($targets eq '0') {$insert = "";}
open my $imp_pairs, '-|', "samtools view -F 2 $insert $bam | awk '\$2 ~ /81|161|97|145|65|129|113|177/'" or die;

my $count = 0;

my $addchrprefix = 1; #if 1, check for "chr" prefix; if not present, add it
if($disablechrprefix == 1) {$addchrprefix = 0;}

my %prevgene = (); #store genes already seen [format: chr_start]

my %prevread = (); #store reads already seen

my %fusepairs = (); #store unique read pairs per fusion

my %genes2exons = (); #store exon coordinates and chromosome for each candidate fusion gene

my %storeids4blast = (); #store read ids of properly paired reads that passed filtration criteria; used to build blast database

my @progress = ('|','/','-','\\');
my $progressItor = 0;

if($VERBOSE == 1) {
    my $currtime = time - $now;
    print "\n   Time                 Current step\n";
    print "---------- --------------------------------------\n";
    printf("[%02d:%02d:%02d]", int($currtime / 3600), int(($currtime % 3600) / 60),int($currtime % 60));
    print " Analyzing improperly paired reads...\n";
}

while(<$imp_pairs>){
    my $line = $_;

    my @tokens = split("\t", $line);


	$count++;
	if($count % 1000 == 0 && $VERBOSE == 1) {
	    print "           $progress[$progressItor] $count reads processed     \r";
	    $progressItor++;
	    if($progressItor > 3) {$progressItor = 0;}
	}
        
	my $id = $tokens[0];

	#skip read ids that have already been processed
	if(!exists($prevread{$id})) {$prevread{$id}=$id;}
	else {
	    next;
	}

	#read1
	my $chr = $tokens[2];
    if($addchrprefix == 1) {$chr = addprefix($chr);}
	my $pos = $tokens[3];

        my $gene = "";
        my $gene_start = "";
        my $gene_end = "";

	#read 2
	my $chrb = $tokens[6];
    if($addchrprefix == 1) {$chrb = addprefix($chrb);}
	if($chrb =~ m/=/) {$chrb = $chr;}
	my $posb = $tokens[7];

        my $geneb = "";
        my $gene_startb = "";
        my $gene_endb = "";

	if(!exists($prevgene{"$chr $pos"})){

#-----------------find closest start coordinate to read1 (and corresponding end)--------------

	    findGene($pos, \@{$exons_start[$chrom{$chr}]}, $gene_start, $gene_end, $gene, $chr);
        
	    $prevgene{"$chr $pos"} = "$gene|$gene_start|$gene_end";#store gene                                                                    
	
	}else{
	    my $tmp = $prevgene{"$chr $pos"};
	    my @vars= split("\\|", $tmp);
	    $gene = $vars[0];
	    $gene_start = $vars[1];
	    $gene_end = $vars[2];
	}

#----------------find closest start coordinate to read2 (and corresponding end)------------

	if(!exists($prevgene{"$chrb $posb"})){

	    findGene($posb, \@{$exons_start[$chrom{$chrb}]}, $gene_startb, $gene_endb, $geneb, $chrb);

	    $prevgene{"$chrb $posb"} = "$geneb|$gene_startb|$gene_endb"; #store gene
	
	   }else{
	       my $tmp = $prevgene{"$chrb $posb"};
	       my @vars = split("\\|", $tmp);
	       $geneb = $vars[0];
	       $gene_startb = $vars[1];
	       $gene_endb = $vars[2];
	   }

#-----------------if candidate translocation------------------------------------------------

	if($gene ne $geneb){ #discordant gene pair...

        
	    next if $gene eq "" || $geneb eq "";

	    my $dist = 0;
	    if($chr eq $chrb){ #calculate distance if same chromosome
		$dist = abs($gene_start - $gene_startb);
	    }

	    my $fusion = "$gene-$geneb";
	 
	    print output "$gene\t$geneb\t$fusion\t$chr\t$chrb\t$gene_start\t$gene_startb\t$dist\t$pos\t$posb\n";
	    
	    #store chromosome and start position of each exon corresponding to candidate fusion gene
	    my $minstart = max(1,min($gene_start, $pos-300));
	    my $minstartb = max(1,min($gene_startb, $posb-300));
        my $maxend = max($gene_end, $pos+300);
        my $maxendb = max($gene_endb, $posb+300);
	    if(!exists($genes2exons{$gene}->{"$chr $minstart"})){ $genes2exons{$gene}->{"$chr $minstart"} = $maxend;}
        elsif($genes2exons{$gene}->{"$chr $minstart"} < $maxend){
            $genes2exons{$gene}->{"$chr $minstart"} = $maxend;
        }
	    if(!exists($genes2exons{$geneb}->{"$chrb $minstartb"})){ $genes2exons{$geneb}->{"$chrb $minstartb"} = $maxendb;}
        elsif($genes2exons{$geneb}->{"$chrb $minstartb"} < $maxendb){
            $genes2exons{$geneb}->{"$chrb $minstartb"} = $maxendb;
        }

        #eliminate redundant fragments (same start for both reads) for depth statistics 

	    if(!exists($fusepairs{"$gene $geneb"}->{"$pos $posb"})){ 
		$fusepairs{"$gene $geneb"}->{"$pos $posb"} = "$pos $posb";
		$storeids4blast{$id} = $id; #store read ids for blast database
	    }        
        
        #if($makeblastdb == 1) {print blast ">$tokens[0] " . ($blastitor++) . "\n$tokens[9]\n";} #print to blast database
        
        }
        
}

close($imp_pairs);
close(output);
if($VERBOSE == 1) {print "\n";}

#-------------print out depth statistics for mapped improper pairs-----------------------


my %genes = (); #store unique genes present in putative translocations

my %pairs = (); #store unique gene pairs present in putative translocations

my %depth = (); #for sorting depth

my $fusetargets = basename($bam);
$fusetargets =~ s/bam/factera.fusiontargets.bed/g;
$fusetargets = "$OUTPUTDIR/$fusetargets";

open(fusiontargets, ">$fusetargets") or die $!;

foreach my $g (keys %fusepairs){
    my $d = (keys %{$fusepairs{$g}});
    $depth{$g} = $d;
}
$count = 0;
foreach my $g (sort {$depth{$b} <=> $depth{$a}} keys %depth){
    $count++;
    next if $depth{$g} < $MINIMPROPERREADS; #skip if less than x discordant reads
    my @gs = split(" ", $g);
    if(!exists($genes{$gs[0]})) {
	$genes{$gs[0]} = $gs[0];
	#print bed coordinates of candidate fusion gene 1
	my $chr = "";
	my $min = 9999999999;
	my $max = 0;
	foreach my $key (keys %{$genes2exons{$gs[0]}}){
	    my @vars = split(" ",$key);
	    $chr = $vars[0];
        if($addchrprefix == 1) {$chr =~ s/chr//g;}
	    my $start = $vars[1];
	    if($start < $min) {$min = $start;}
	    my $end = $genes2exons{$gs[0]}->{$key};
	    if($end > $max){$max = $end;}
	}
	print fusiontargets "$chr\t$min\t$max\t$g\n";
    }
    if(!exists($genes{$gs[1]})) {
	$genes{$gs[1]} = $gs[1];
	#print bed coordinates of candidate fusion gene 2 //TODO: extend to maximum read start (after last/first exon)
	my $chr= "";
	my $min= 9999999999;
	my $max= 0;
        foreach my $key (keys %{$genes2exons{$gs[1]}}){
            my @vars = split(" ",$key);
            $chr = $vars[0];
            if($addchrprefix == 1) {$chr =~ s/chr//g;}
            my $start = $vars[1];
            if($start <$min) {$min = $start;}
            my $end = $genes2exons{$gs[1]}->{$key};
            if($end > $max){$max = $end;}
        }
        print fusiontargets "$chr\t$min\t$max\t$g\n";
    }
    print output2 "$g\t$depth{$g}\n";
    $pairs{$gs[0]}->{$gs[1]} = $depth{$g}; #store fusion partner of gene 1
    $pairs{$gs[1]}->{$gs[0]} = $depth{$g}; #store fusion partner of gene 2
}

close(fusiontargets);
close(output2);

#===============================================================================================
#------------infer break-point by soft-clipping, and increment depth of fusion pairs accordingly

$count = 0; #reset counter for progress meter

#open my $sclipped, '-|', "samtools view -F 12 -L $fusetargets $bam | awk '\$2 ~ /99|147|83|163|67|131|115|179/' | awk '\$6 ~ /S/'" or die;

open my $sclipped, '-|', "samtools view -F 12 -L $fusetargets $bam | awk '\$2 ~ /99|147|83|163|67|131|115|179|81|161|97|145|65|129|113|177/' | awk '\$6 ~ /S/'" or die;

my %breakpoints = (); #store potential breakpoints [key1=gene; key2=position; value=count]

my %screads = (); #store soft-clipped reads meeting gene and length criteria

my %consensus = (); #store consensus sequence for each breakpoint cluster

$addchrprefix = 1; #if 1, check for "chr" prefix; if not present, add it

if($VERBOSE == 1) {
    my $currtime = time - $now;
    printf("[%02d:%02d:%02d]", int($currtime / 3600), int(($currtime % 3600) / 60),int($currtime % 60));
    print " Analyzing soft-clipped reads...\n";
}

while(<$sclipped>){
    my $line = $_;

    my @tokens = split("\t", $line);

    my $cigar = $tokens[5];
   
    #next if !($cigar =~ m/S/);
    next if ($cigar =~ m/D/);
    next if($cigar =~ m/I/);
 
    $count++;
    if($count % 100 == 0 && $VERBOSE == 1) {
        print "           $progress[$progressItor] $count reads processed     \r";
        $progressItor++;
        if($progressItor > 3) {$progressItor = 0;}
    }

	my $tmp = $cigar;
	$tmp =~ s/S//g;
	next if length($tmp) != length($cigar) - 1;

	my @parsed = parsecigar($cigar);
    my @parsed2 = parsecigar($parsed[2]);

	my $matched = 0; #cigar string 'M'
	my $skipped = 0; #cigar string 'S'

	#if forward orientation, properly paired, and with match followed by clip, go to next..
	if(($tokens[1] == 99) || ($tokens[1] == 163)){
	    next if $parsed[1] eq "M"; 
	}
	#if reverse orientation, properly paired, and with clip followed by match, go to next..                                                      
	if(($tokens[1] == 147) || ($tokens[1] == 83)){
            next if $parsed[1] eq "S";
        }

	#only keep clipped regions at least clipsize (default 16 bases; only 1 in 4.3B by random chance)
	if($parsed[1] eq "S") {$skipped = $parsed[0]; $matched = $parsed2[0]; next if $parsed[0] < $clipsize;}
	elsif($parsed2[1] eq "S"){$skipped = $parsed2[0]; $matched = $parsed[0]; next if $parsed2[0] < $clipsize;}

	my $chr = $tokens[2];
    if($addchrprefix == 1) {$chr = addprefix($chr);}
    my $pos = $tokens[3];
	my $gene = "";
	my $gene_start = 0;
	my $gene_end = 0;

        #find gene

	findGene($pos, \@{$exons_start[$chrom{$chr}]}, $gene_start, $gene_end, $gene, $chr);

	next if !exists($genes{$gene});
        
    if($makeblastdb == 1) {
        #so makeblastdb does not complain about >40% Ns
        my $Ns = $tokens[9];
        $Ns =~ s/N//g;
        if(length($Ns)/length($tokens[9]) >= .75) {
            print blast ">$tokens[0]\_" . ($blastitor++) . "\_SC\n$tokens[9]\n";
        }
    } #print to blast database

	my $bp = $pos;
	if($parsed[1] eq "M") {
	    $bp += $matched - 1;

	    my $prevLen = @{$screads{$gene}->{$bp}}[0];
	    my $clipped = substr($tokens[9],$matched);
	    my $notclipped = substr($tokens[9],0,$matched);
	    my $middle = length($tokens[9])/2;
        
        
	    if(abs(length($clipped) - $middle) <= abs($prevLen - $middle)){
            my @temp = (length($clipped), "$notclipped $clipped", "NC", $tokens[0]);
            $screads{$gene}->{$bp} = \@temp;
            my $adjustlen = 500;            
            #store consensus sequence for soft-clipped portion of a given breakpoint
            my @seqarray = getSeqArray(\%consensus, $gene, $bp, 0, $adjustlen - $clipsize, 0, length($clipped), \%screads, $clipped);
            $consensus{$gene}->{$bp} = \@seqarray;
            
	    }
	}else{
	    my $prevLen = @{$screads{$gene}->{$bp}}[0];
        my $clipped = substr($tokens[9],0,length($tokens[9])-$matched);
	    my $notclipped = substr($tokens[9],length($tokens[9])-$matched);
	    my $middle = length($tokens[9])/2;
	    if(abs(length($clipped) - $middle) <= abs($prevLen - $middle)){
            my @temp = (length($clipped), "$clipped $notclipped", "CN", $tokens[0]);
            $screads{$gene}->{$bp} = \@temp;
         
            #store consensus sequence for soft-clipped portion of a given breakpoint
	    my $adjustlen = 500;
            my @seqarray = getSeqArray(\%consensus, $gene, $bp, $adjustlen - $clipsize - length($clipped), $adjustlen - $clipsize, 0, length($clipped), \%screads, $clipped);
            $consensus{$gene}->{$bp} = \@seqarray;
	    }

	}

	$breakpoints{$gene}->{"$chr:$bp"}++;
    
}

if($VERBOSE == 1) {print "\n";}

#------------------------------------------------------------------------------
#finish writing reads to disk, then make blast database if not already made
if($makeblastdb == 1) {makeBLASTdb();}

#------------------------------------------------------------------------------
#sort soft-clipped reads by breakpoint, identify "best" breakpoints and analyze

my %usedpair = (); #store gene pairs to output only unique pairs

my %bestfusions = (); #only keep fusions with best support for a given gene pair and pair of coordinates (+/-5 bp pad)

my $fusions = basename($bam);
$fusions =~ s/bam/factera.fusions.txt/g;
$fusions = "$OUTPUTDIR/$fusions";

open(output3, ">$fusions") or die $!;

my $fusionseqs = basename($bam);
$fusionseqs =~ s/bam/factera.fusionseqs.fa/g;
$fusionseqs = "$OUTPUTDIR/$fusionseqs";

open(output4, ">$fusionseqs") or die $!;

my $fusions_bed = basename($bam);
$fusions_bed =~ s/bam/factera.fusions.bed/g;
$fusions_bed = "$OUTPUTDIR/$fusions_bed";

open(output5, ">$fusions_bed") or die $!;

if($VERBOSE == 1) {
    my $currtime = time - $now;
    printf("[%02d:%02d:%02d]", int($currtime / 3600), int(($currtime % 3600) / 60),int($currtime % 60));
    print " Validating candidate fusions...\n";
}

#==================================================================================================
my $total_comp = 0; #count total events to be compared (for progressbar)
if($VERBOSE == 10) {
    foreach my $gene (keys %breakpoints){
        my $itor = 0;
        foreach my $bp (sort {$breakpoints{$gene}->{$b} <=> $breakpoints{$gene}->{$a}} keys %{$breakpoints{$gene}}){
            $itor++;
            last if $itor > $MAXBPS2EXAMINE; #only examine the x most abundant putative breakpoints for gene 1
            my $count = $breakpoints{$gene}->{$bp};
            next if $count < $MINBPSUPPORT; #skip if <x reads support
             my $itor2 = 0;
            foreach my $gene2 (keys %{$pairs{$gene}}) {
                foreach my $bp2 (sort {$breakpoints{$gene2}->{$b} <=> $breakpoints{$gene2}->{$a}} keys %{$breakpoints{$gene2}}){
                    $itor2++;
                    last if $itor2 > $MAXBPS2EXAMINE; #examine x most abundant putative breakpoints for gene 2
                    $total_comp++;
                }
            }
        }
    }
}
#==================================================================================================

my $progress_itor = 0;
$progressItor = 0;
my $prev_prog = -1;
my $fusioncons = "";

foreach my $gene (keys %breakpoints){
    
    my %rankedlist = ();
    #Rank all breakpoints for partner genes wrt read depth
    foreach my $gene2 (keys %{$pairs{$gene}}) {
        
        foreach my $bp2 (sort {$breakpoints{$gene2}->{$b} <=> $breakpoints{$gene2}->{$a}} keys %{$breakpoints{$gene2}}){
            my $depth = $breakpoints{$gene2}->{$bp2};
            $rankedlist{$depth}->{$gene2}->{$bp2} = $bp2;
        }
    }
    
    
    my $itor = 0;
    foreach my $bp (sort {$breakpoints{$gene}->{$b} <=> $breakpoints{$gene}->{$a}} keys %{$breakpoints{$gene}}){
        $itor++;
        last if $itor > $MAXBPS2EXAMINE; #only examine the x most abundant putative breakpoints for gene 1
        my $count = $breakpoints{$gene}->{$bp};
        next if $count < $MINBPSUPPORT; #skip if <x reads support

        my $bp_ = $bp;
        $bp_ =~ s/.*://g;
        my $chr = $bp;
        $chr =~ s/:.*//g;
        
        my $id = @{$screads{$gene}->{$bp_}}[3];
     
        my $clipOrder = @{$screads{$gene}->{$bp_}}[2];
        my $read1 = @{$screads{$gene}->{$bp_}}[1];
        my $read1b = $read1;
        my $readLen = length($read1) - 1;
        if($READLENBIT == 0) {$READLEN = $readLen; $READLENBIT = 1;}
       
        my $read1cut = @{$screads{$gene}->{$bp_}}[0];
        if($clipOrder eq "NC"){ #excise non-clipped read using order bit     
            $read1 = substr($read1, 0, (length($read1)-$read1cut)-1);
            my $read1btemp = substr($read1b, length($read1b)-$read1cut);
            
            $read1b = getConsensus(\%consensus, $gene, $bp_);
	    $read1b = substr($read1b, 0, length($read1btemp));
#	    $read1b = $read1btemp;
        }else{
            $read1 = substr($read1, ($read1cut + 1));
            my $read1btemp = substr($read1b, 0, $read1cut);
            
            $read1b = getConsensus(\%consensus, $gene, $bp_);
	    $read1b = substr($read1b, 1 + length($read1b) - length($read1btemp));
#	    $read1b = $read1btemp;
        }
        my %kmers1 = getKmers($read1, $k);
        my %kmers1rc = getKmers(doRevComp($read1), $k);
        my %kmers1b = getKmers($read1b, $k);
        my %kmers1rcb = getKmers(doRevComp($read1b), $k);
  
        
        my $itor2 = 0;
       
        foreach my $depth (sort {$b <=> $a} keys %rankedlist) {
        foreach my $gene2 (keys %{$rankedlist{$depth}}) {
            foreach my $bp2 (keys %{$rankedlist{$depth}->{$gene2}}){
                $itor2++;
                last if $itor2 > $MAXBPS2EXAMINE; #examine x most abundant putative breakpoints paired with gene 1
                my $count2 = $breakpoints{$gene2}->{$bp2};
                next if $count2 < $MINBPSUPPORT; #skip if <x reads support
                
                
                if($VERBOSE == 10){
                    $progress_itor++;
                    my $prog = int(100*$progress_itor/$total_comp);
                    if($prog % 4 == 0 && $prog != $prev_prog){
                        print "*";
                        $prev_prog = $prog;
                    }
                }
                $progress_itor++;
                if($progress_itor % 100 == 0 && $VERBOSE == 1) {
                    print "           $progress[$progressItor] $progress_itor candidates processed     \r";
                    $progressItor++;
                    if($progressItor > 3) {$progressItor = 0;}
                }
                
                my $bp2_ = $bp2;
                $bp2_ =~ s/.*://g;
                my $chr2 = $bp2;
                $chr2 =~ s/:.*//g;
                my $id2 = @{$screads{$gene2}->{$bp2_}}[3];
   
                my $clipOrder2 = @{$screads{$gene2}->{$bp2_}}[2];
                my $read2 = @{$screads{$gene2}->{$bp2_}}[1];
                my $read2b = $read2;

                my $read2cut = @{$screads{$gene2}->{$bp2_}}[0];
                if($clipOrder2 eq "NC"){ #excise clipped read using order bit
                    my $read2temp = substr($read2,length($read2)-$read2cut);
                    $read2b = substr($read2b, 0, (length($read2b)-$read2cut)-1);
                    
                    $read2 = getConsensus(\%consensus, $gene2, $bp2_);
		    $read2 = substr($read2, 0, length($read2temp));
                }else{
                    my $read2temp = substr($read2, 0, $read2cut);
                    $read2b = substr($read2b,($read2cut+1));
                    
                    $read2 = getConsensus(\%consensus, $gene2, $bp2_);
		    $read2 = substr($read2, 1 + length($read2) - length($read2temp));
                }
                my $cmpthresh = 0.5;
                #compare clipped region of read2 with read1 non-clipped
                my $threshold = min(25,$cmpthresh*(min(length($read2), (length($read1))) - $k));
                my $thresholdb = min(25,$cmpthresh*(min(length($read2b), (length($read1b))) - $k));
                my $doOffset = 1; #if 1, allow offset adjustment

                
                my @same = compKmers(\%kmers1, \%kmers1rc, $read2, $k, $threshold, $read1, $doOffset, $clipOrder, $clipOrder2);
                my $clipOrderb = "NC";
                my $clipOrder2b = "NC";
                if($clipOrder eq "NC") {$clipOrderb = "CN";}
                if($clipOrder2 eq "NC") {$clipOrder2b = "CN";}
                my @sameb = compKmers(\%kmers1b, \%kmers1rcb, $read2b, $k, $thresholdb, $read1b, $doOffset, $clipOrderb, $clipOrder2b);
                
                if($same[0] != 9999999 && $same[0] != $sameb[0]) {$same[0] = $sameb[0] = 9999999;}
                
                #if clipped region of read 2 is sufficiently similar to non-clipped region of read 1...
                if($same[0] != 9999999 && $sameb[0] != 9999999){
                    
                    next if exists($usedpair{"$gene\_$gene2\_$bp\_$bp2\_$clipOrder\_$clipOrder2"}) || exists($usedpair{"$gene2\_$gene\_$bp2\_$bp\_$clipOrder2\_$clipOrder"});
                    next if ($clipOrder eq $clipOrder2) && $same[1] eq 'F'; #'F' = forward read
                    next if($clipOrder ne $clipOrder2) && $same[1] eq 'RC'; #'RC' = reverse complement read 
                    $usedpair{"$gene\_$gene2\_$bp\_$bp2\_$clipOrder\_$clipOrder2"} = 1; #store unique gene pair

                    #assign strand orientation to fusion
                    my $orientation1 = "2-";
                    if($clipOrder eq "NC") {$orientation1 = "1+";}
                    my $orientation2= "2-";
                    if($clipOrder2 eq "CN") {$orientation2 = "1+";}
                    if(($clipOrder ne $clipOrder2) && $clipOrder eq "CN") {$orientation1 = "2+"; $orientation2 = "1+";}
                    if(($clipOrder ne $clipOrder2) && $clipOrder eq "NC") {$orientation2 = "2+"; $orientation1 = "1+";}
                    
                    
                    #do bp correction by comparison to reference=================

                    my $offsetbp1 = 0;
                    my $offsetbp2 = $same[0];
                    
                    if($same[0] != 0){
                        my @tmp = bp_correction(@{$screads{$gene2}->{$bp2_}}[1], \@same, $offsetbp2, $bp2_, $chr2);
                        $offsetbp1 = $tmp[0];
                        $offsetbp2 = $tmp[1];
                    }
                    
                    #=============================================================

                    #Find and store insert (non-templated or microhomologous) sequence (if any)
                    if($offsetbp2!=0){
                        if($offsetbp2>0 && $clipOrder2 eq "NC") {$fusioncons = substr($read2, 0, $offsetbp2);}
                        elsif($offsetbp2<0 && $clipOrder eq "NC") {$fusioncons = doRevComp(substr($read1b, 0, -$offsetbp2));}
                    }
        
                    
     
                    #return and print adjusted fragment 2 breakpoint
                    
                    my $fusionSeqFa = getFusionSeq($clipOrder, $clipOrder2, $bp_, $bp2_, $bp, $bp2, $offsetbp2, $buffer, $offsetbp1);
                    
                    my @vals = split("\n", $fusionSeqFa);
                    my $fusionSeq = $vals[1];
                    my $part1 = substr($fusionSeq, 0, $buffer);
                    my $part2 = substr($fusionSeq, $buffer, $buffer);
                    #print "$part1\t$part2\n";
                    #print "$fusioncons\t$bp_\t$bp2_$orientation1\t$orientation2\n";
                    #add non-templated sequence into estimated fusion sequence====
                    my $inserttype = "-";
                    if(length($fusioncons) > 0){
                        my $nt_len = length($fusioncons);
                        if($orientation1 =~ m/1/) {
                            $fusioncons = doRevComp($fusioncons);
                            if(substr($part2, 0, length($fusioncons)) eq $fusioncons) {
                                #$inserttype = "Microhomology";
                                $fusioncons = "";
                            }
                            else {
                                
                                $inserttype = "Non-templated";
                                substr($part2, 0, length($fusioncons)) = $fusioncons;
                                
                                if($orientation2 =~ m/\+/){
                                    $bp2_ += $nt_len;
                                }else {$bp2_ -= $nt_len;}
                                
                            }
                            
                        } else {
                            if(substr($part1, length($part1)-length($fusioncons), length($fusioncons)) eq $fusioncons){
                                #$inserttype = "Microhomology";
                                $fusioncons = "";
                            }else {
                                $inserttype = "Non-templated";
                                substr($part1, length($part1)-length($fusioncons), length($fusioncons)) = $fusioncons;
                                
                                if($orientation1 =~ m/\+/){
                                    $bp2_ -= $nt_len;
                                }else {$bp2_ += $nt_len;}
                            }
                        }
                    }
                    
                    #print "$fusioncons\t$bp_\t$bp2_\n";
                    
                    if($fusioncons eq "") {$fusioncons = "-"; $fusioncons_bp = 0;}
                    else {$fusioncons = "[$fusioncons]";}
                    $nonref_insert = "";
                    $fusionSeqFa = $vals[0] . "\n" . $part1 . "" . $part2;
                    #=============================================================
                    
                    my $data = "$gene\t$gene2\t$chr:" . ($bp_ + $offsetbp1) . "\t$chr2:" . ($bp2_ + $offsetbp2) . "\t$count\t$count2\t$same[0]\t$orientation1 $orientation2\t$clipOrder\t$clipOrder2\t";
                    
                    #do BLAST search
                    my $datatmp = doBLAST($fusionSeqFa, $gene, $gene2, $buffer, ($MINBLASTFRAC * $READLEN));
                    $data = $data . "$datatmp";
                    
                    #print close-up of breakpoint
                    
                    my @vals = split("\n", $fusionSeqFa);
                    my $fusionSeq = $vals[1];
                    $part1 = substr($fusionSeq, $buffer - $pad, $pad);
                    $part2 = substr($fusionSeq, $buffer, $pad);
                    
                    #add non-templated sequence into estimated fusion sequence

                    
                    $data = $data . $part1 . " <> " . $part2 . " $fusioncons";#\_$inserttype";
                    $fusioncons = "";
                    #print "$data\n";
                    my @vars = split("\t", $datatmp);
                    my $bp_depth = $vars[0]; #retrieve breakpoint depth
                    next if $bp_depth < $MINSPANNINGREADS; #don't save putative translocations with too few reads spanning breakpoint
                    my @tmp = ($data, $fusionSeqFa);
                    my $bp1breaks = $tmp[4];
                    my $bp2breaks = $tmp[5];
                    my $code = "$chr:" . ($bp_ + $offsetbp1) . " $chr2:" . ($bp2_ + $offsetbp2);
                    #only keep best fusion with same coordinates and depth (highest read 1 and read 2 bp support)
                    if(!exists($bestfusions{$bp_depth}->{$code})) {$bestfusions{$bp_depth}->{$code} = \@tmp;
                    }else{
                        my @tmp2 = $bestfusions{$bp_depth}->{$code};
                        if($bp1breaks > $tmp2[4] && $bp2breaks > $tmp2[5]){
                            $bestfusions{$bp_depth}->{$code} = \@tmp;
                        }
                    }
                }
            }
        }
    }
}
    
}

#if(($progress_itor == $total_comp || $progress_itor == 0) && $VERBOSE == 1) {print "\n";}

#print header
if(keys %bestfusions > 0){
    print output3 "Est_Type\tRegion1\tRegion2\tBreak1\tBreak2\tBreak_support1\tBreak_support2\tBreak_offset\tOrientation\tOrder1\tOrder2\tBreak_depth\tProper_pair_support\tUnmapped_support\tImproper_pair_support\tPaired_end_depth\tTotal_depth\tFusion_seq\tNon-templated_seq\n";
    if($VERBOSE == 1) {print "\n-------------------------------------------------\n\nFACTERA results:\nEst_Type\tRegion1\tRegion2\tBreak1\tBreak2\tBreak_support1\tBreak_support2\tBreak_offset\tOrientation\tOrder1\tOrder2\tBreak_depth\tProper_pair_support\tUnmapped_support\tImproper_pair_support\tPaired_end_depth\tTotal_depth\tFusion_seq\tNon-templated_seq\n";}
}

#===========================================================================================================       
#print final fusion predictions, ordered by decreasing breakpoint depth (remove redundancies with 20bp window)
my %bestfusions_ = (); #store best fusion events and use to remove redundancies (+/- 20bp padding)

my %normdepth = getNormDepth(\%bestfusions, $buffer); #total depth of genomic regions flanking candidate fusion(s) 

my $transitor = 0; #count translocations/fusions
foreach my $bp_depth (sort {$b <=> $a} keys %bestfusions){
    foreach my $bps (keys %{$bestfusions{$bp_depth}}){
        my $orig_bps = $bps;
        my @tmp = @{$bestfusions{$bp_depth}->{$bps}};
        my $data = $tmp[0];
        my @tokens = split(" ", $bps);
        my $bp1 = $tokens[0];
        my $bp2 = $tokens[1];
        my $chr1 = $bp1;
        $chr1 =~ s/:.*//g;
        $bp1 =~ s/.*://g;
        my $chr2 = $bp2;
        $chr2 =~ s/:.*//g;
        $bp2 =~ s/.*://g;
        my @vars = split("\t", $data);
        my $orient1 = $vars[7];
        
        #Estimate rearrangement type
        my $first = substr($orient1,0,1);
        my $second = substr($orient1,3,1);
        my $first_polarity = substr($orient1,1,1);
        my $sec_polarity = substr($orient1,4,1);
        
        my $SR_type = " - ";
        if($chr1 ne $chr2) {$SR_type = "TRA";}
        elsif($first_polarity ne $sec_polarity) {$SR_type = "INV";}
        elsif(($bp1<$bp2 && $first == 1) || ($bp1>$bp2 && $first == 2)){$SR_type = "DEL";}
        #elsif(abs(1+$bp1-$bp2)>1000000){$SR_type = "TRA";}
        
        #Switch order for instances where region 2 is first (5' of fusion)
        if($first == 2){
            my $tmp = $vars[0];
            $vars[0] = $vars[1];
            $vars[1] = $tmp;
            $tmp = $vars[2];
            $vars[2] = $vars[3];
            $vars[3] = $tmp;
            $tmp = $vars[4];
            $vars[4] = $vars[5];
            $vars[5] = $tmp;
            $tmp = $vars[8];
            $vars[8] = $vars[9];
            $vars[9] = $tmp;
            $vars[7] =~ s/1/0/g;
            $vars[7] =~ s/2/1/g;
            $vars[7] =~ s/0/2/g;
            $orient1 = $vars[7];
            $tmp = $bp1;
            $bp1 = $bp2;
            $bp2 = $tmp;
            $tmp = $chr1;
            $chr1 = $chr2;
            $chr2 = $tmp;
            $bps = "$chr1:$bp1 $chr2:$bp2";
        }
        
        #remove redundancies within 20bp window
        foreach my $bp_depth_ (sort {$b <=> $a} keys %bestfusions){
            foreach my $bps_ (keys %{$bestfusions{$bp_depth_}}){
                my @tokens2 = split(" ", $bps_);
                my $bp1_ = $tokens2[0];
                my $bp2_ = $tokens2[1];
                my $chr1_ = $bp1_;
                $chr1_ =~ s/:.*//g;
                $bp1_ =~ s/.*://g;
                my $chr2_ = $bp2_;
                $chr2_ =~ s/:.*//g;
                $bp2_ =~ s/.*://g;
                my @tmp_ = @{$bestfusions{$bp_depth_}->{$bps_}};
                my $data_ = $tmp_[0];
                my @vars_ = split("\t", $data_);
                
                my $orient2 = $vars_[7];
                my $first = substr($orient2,0,1);
                
                #Switch order for instances where region 2 is first (5' of fusion)
                if($first == 2){
                    $vars_[7] =~ s/1/0/g;
                    $vars_[7] =~ s/2/1/g;
                    $vars_[7] =~ s/0/2/g;
                    $orient2 = $vars_[7];
                    my $tmp = $bp1_;
                    $bp1_ = $bp2_;
                    $bp2_ = $tmp;
                    $tmp = $chr1_;
                    $chr1_ = $chr2_;
                    $chr2_ = $tmp;
                    $bps_ = "$chr1_:$bp1_ $chr2_:$bp2_";
                }
                next if $bps eq $bps_;
                next if ($chr1 ne $chr1_) || ($chr2 ne $chr2_);
                next if $orient1 ne $orient2;

                if(($bp1_ <= $bp1 + 20) && ($bp1_ >= $bp1 - 20) && ($bp2_ <= $bp2 + 20) && ($bp2_ >= $bp2 - 20)){
                    $bestfusions_{$bps_} = $bps_;
                }
            }
        }
        
        if(!exists($bestfusions_{$bps})) {
            $bestfusions_{$bps} = $bps;
            my $seq = $vars[(@vars-1)];
            my @tokens = split(" ", $seq);
            $seq = join("\t", @tokens);
            $vars[(@vars-1)] = $normdepth{$orig_bps}; #retrieve genomic depth
            $data = join("\t", @vars);
            $data = $SR_type . "\t" . $data . "\t". $seq;
            print output3 "$data\n";
            if($VERBOSE == 1) {print "$data\n";}
            my $fusionSeq = $tmp[1];
            print output4 "$fusionSeq\n";
            my @t1 = split(":", $vars[2]);
            my @t2 = split(":", $vars[3]);
            print output5 "$t1[0]\t$t1[1]\t$t1[1]\t$vars[0]\_$vars[1]\_$vars[2]\_$vars[3]\n";
            print output5 "$t2[0]\t$t2[1]\t$t2[1]\t$vars[0]\_$vars[1]\_$vars[2]\_$vars[3]\n";
            $transitor++; #count translocations
        }
    }
}

if($VERBOSE == 1) {
    if($transitor > 0) {
        my $ins = "fusion";
        if($transitor>1) {$ins = "fusions";}
        print "$transitor $ins found!\n";
    }
    else {
        print "\nNo fusions found\n";
    }
}

close(output3);
close(output4);
close(output5);

$now = time - $now;
# Print runtime #

if($TIME == 1) {printf("Total running time: %02d:%02d:%02d\n", int($now / 3600), int(($now % 3600) / 60),int($now % 60));}




#===========================================================================================================      
#FUNCTIONS................

#===========================================================================================================
#retrieve depth of region surrounding fusion(s)

sub getNormDepth{
    
    my $buffer = $_[1];

    my %coors = (); #store surrounding depth for each translocation

    foreach my $bp_depth ( keys %{$_[0]}){
        foreach my $bps (keys %{$_[0]{$bp_depth}}){
    
            my $data = $_[0]{$bp_depth}->{$bps};
            my @tokens = split(" ", $bps);
            my $bp1 = $tokens[0];
            my $bp2 = $tokens[1];
            my $chr1 = $bp1;
            if($addchrprefix == 1) {$chr1 =~ s/chr//g;}
            $chr1 =~ s/:.*//g;
            $bp1 =~ s/.*://g;
            my $chr2 = $bp2;
            $chr2 =~ s/:.*//g;
            if($addchrprefix == 1) {$chr2 =~ s/chr//g;}
            $bp2 =~ s/.*://g;
            my @vars = split("\t", $data);
            my $orient1 = $vars[7];
            my @cut = split(" ", $orient1);
            if($cut[0] =~ m/-/){
                $bp1 -= $buffer;
            }
            $coors{$chr1}->{$bp1} = $tokens[0]; #chr1:bp1
            if($cut[1] =~ m/-/){
                $bp2 -= $buffer;
            }
            $coors{$chr2}->{$bp2} = $tokens[1]; #chr2:bp2
   
            my $stat = Statistics::Descriptive::Full->new();
            
            my $insert = "-l $targets";
            if($targets eq '0') {$insert = "";}

            my $pid = open3(gensym, \*PH, ">&NULL", "samtools mpileup -ABr $chr1:$bp1-" . ($bp1+$buffer) . " -Q20 $insert -d 10000000 $bam");

            my @d1 = (); #store bp depth for gene 1
            
            while(<PH>){
                my $line = $_;
                chomp($line);
                my @vars = split("\t", $line);
                my $depth = $vars[3];
                push(@d1, $depth);
            }
            
            close(PH);
            
            $stat->add_data(@d1);
            my $med_depth1 = 0;
            $med_depth1 = $stat->median();

            
            my $pid = open3(gensym, \*PH, ">&NULL", "samtools mpileup -ABr $chr2:$bp2-" . ($bp2+$buffer) . " -Q20 $insert -d 10000000 $bam");

            my @d2 = (); #store bp depth for gene 2
            
            while(<PH>){
                my $line = $_;
                chomp($line);
                my @vars = split("\t", $line);
                my $depth = $vars[3];
                push(@d2, $depth);
            }
            
            close(PH);
            
            $stat = Statistics::Descriptive::Full->new();
            $stat->add_data(@d2);
            my $med_depth2 = 0;
            $med_depth2 = $stat->median();
            
            my $final_depth = 0;
            
            #if(@d1 > 0 && @d2 > 0){
                # $final_depth = ($med_depth1 + $med_depth2)/2;
                #}else {
                $final_depth = max($med_depth1, $med_depth2);
                #}

            $coors{$bps} = $final_depth;
        }
        
    }
    return %coors;

}

#===========================================================================================================       
#BLAST all reads that map to genes involved in fusion region (properly paired, improperly paired, and not paired)

sub doBLAST{
       
    my $gene = $_[1];
    my $gene2 = $_[2];
    my $buffer = $_[3];
    my $min_len = $_[4];
    
    #write fusion gene to temp file for blast search
    my $query = basename($bam);
    $query =~ s/bam/factera.blastquery.fa/g;
    $query  = "$OUTPUTDIR/$query";
    open(bquery, ">$query") or die $!;
    print bquery "$_[0]";
    close(bquery);
    
    my $bdbname = basename($bam);
    $bdbname =~ s/bam/factera.blastreads.fa/g;
    $bdbname = "$OUTPUTDIR/$bdbname";
    
    #TODO: detect blast version: blastx: -max_target_seqs; blastn "-num_alignments 9999999 -num_descriptions 9999999"
    #open my $blastout, '-|', "blastn -task 'megablast' -query $query -db $bdbname -outfmt 6 -num_threads $BLASTTHREADS -num_alignments 9999999 -num_descriptions 9999999" or die;
    open my $blastout, '-|', "blastn -task 'megablast' -query $query -db $bdbname -outfmt 6 -num_threads $BLASTTHREADS -max_target_seqs 9999999" or die;
    
    my $bp_depth = 0; #breakpoint depth
    my $bpd_SC = 0; #count soft-clipped properly paired read support
    my $bpd_UM = 0; #count unpaired read support
    my $bpd_IP = 0; #count soft-clipped improperly paired read support
    
    my $proper_pair_depth = 0; #extra support due to properly paired reads that flank breakpoint
    my %properpairs = (); #store proper pairs

    my %storeused = (); #store read ids to ensure each fragment only counted once
    
    my $percentsim = 100 * $minsimilarity;
    
    #iterate through BLAST results and count reads overlapping breakpoint
    while(<$blastout>){
        my $line = $_;
        chomp($line);
        my @tokens = split("\t", $line);
        my $id = $tokens[1];
        $id =~ s/\_.*//g;
	my $readtype = $tokens[1];
	$readtype =~ s/.*\_//g; #store read type code (SC=soft-clipped, UM=unmapped, IP=improperly paired)

        my $percent = $tokens[2];
        my $align_len = $tokens[3];
        my $left = $tokens[6];
        my $right = $tokens[7];

        if($percent >= $percentsim && $align_len >= $min_len){

	    if (exists($storeused{$id})){ #skip if id already used
		$properpairs{$id}->{$left} = $right; #group proper pairs
		next;
	    }

            if(min($left,$right) < ($buffer - 15) && max($left,$right) > ($buffer + 15)){
                $bp_depth++;
                if($readtype eq "SC") {$bpd_SC++;}
                if($readtype eq "UM") {$bpd_UM++;}
                if($readtype eq "IP") {$bpd_IP++;}
                $storeused{$id} = $id; #store id
            }
            $properpairs{$id}->{$left} = $right; #group proper pairs            
        }
    }
  #=====================================  
    #count proper pairs that flank breakpoint
    foreach my $id (keys %properpairs){
        my $size = (keys %{$properpairs{$id}});
        next if $size != 2;

        my @coors = ();
        foreach my $left (keys %{$properpairs{$id}}){
            push(@coors, $left);
            push(@coors, $properpairs{$id}->{$left});
        }
        if((min(@coors) <  ($buffer - 85)) && (max(@coors) >  ($buffer + 85))){
            $proper_pair_depth++;
        }
    }
 #=====================================   
      
    #print break-point depth
    my $data = "$bp_depth\t$bpd_SC\t$bpd_UM\t$bpd_IP\t"; 
    
    #print proper pair depth
    $data = $data . "$proper_pair_depth\t"; 

    return $data;
}

#=========================================================================================================== 
sub makeBLASTdb{
    my $currtime = time - $now;
    printf("[%02d:%02d:%02d]", int($currtime / 3600), int(($currtime % 3600) / 60),int($currtime % 60));
    print" Extracting reads for blast database...\n";
    
    #print all unmapped reads to blast database fasta file 

    open my $blastreads, '-|', "samtools view -F 2 -L $fusetargets $bam" or die; 
    
    $count = 0;
    while(<$blastreads>){
        my $line = $_;
        chomp($line);
        my @tokens = split("\t", $line);        
        if ($storeids4blast{$tokens[0]}){
            $count++;
            if($count % 1000 == 0 && $VERBOSE == 1) {
                print "           $progress[$progressItor] $count reads written     \r";
                $progressItor++;
                if($progressItor > 3) {$progressItor = 0;}
            }
             #so makeblastdb does not complain about >40% Ns
            my $Ns = $tokens[9];
            $Ns =~ s/N//g;
            if(length($Ns)/length($tokens[9]) < .75) {next;}
            print blast ">$tokens[0]\_" . ($blastitor++) . "\_IP\n$tokens[9]\n";
        }elsif(($tokens[1] == 73) || ($tokens[1] == 133) || ($tokens[1] == 89) || ($tokens[1] == 121) || ($tokens[1] == 165) || ($tokens[1] == 181) || ($tokens[1] == 101) || ($tokens[1] == 117) || ($tokens[1] == 153) || ($tokens[1] == 185) || ($tokens[1] == 69) || ($tokens[1] == 137) || ($tokens[1] == 77) || ($tokens[1] == 141)){
            $count++;
            if($count % 1000 == 0 && $VERBOSE == 1) {
                print "           $progress[$progressItor] $count reads written\r";
                $progressItor++;
                if($progressItor > 3) {$progressItor = 0;}
            }
            #so makeblastdb does not complain about >40% Ns
            my $Ns = $tokens[9];
            $Ns =~ s/N//g;
            if(length($Ns)/length($tokens[9]) < .75) {next;}
            print blast ">$tokens[0]\_" . ($blastitor++) . "\_UM\n$tokens[9]\n";
        }
    }
    print "\n";
    close($blastreads);
    
    close(blast);
    
    my $bdbname = basename($bam);
    $bdbname =~ s/bam/factera.blastreads.fa/g;
    $bdbname = "$OUTPUTDIR/$bdbname";
    
    $currtime = time - $now;
    printf("[%02d:%02d:%02d]", int($currtime / 3600), int(($currtime % 3600) / 60),int($currtime % 60));
    print " Creating blast database...";
    
    system("makeblastdb -in $bdbname -dbtype 'nucl' >/dev/null");
    
    if($VERBOSE == 1){
        print "\n           - Done\n";
    }
    
    }



#===========================================================================================================       
#calculate fusion coordinates and retrieve genomic sequences surrounding breakpoint using user-definable pad length (default 200)

sub getFusionSeq{
   
    my $g1start = 0;
    my $g1end = 0;
    my $g2start = 0;
    my $g2end = 0;
    my $buffer = $_[7]; #size of genomic padding around breakpoint
    my $clipOrder = $_[0]; #NC or CN
    my $clipOrder2 = $_[1]; #NC or CN
    my $bp_ = $_[2]; #breakpoint fragment 1
    my $bp2_ = $_[3]; #breakpoint fragment 2
    my $bp = $_[4]; #chr:breakpoint fragment 1
    my $bp2 = $_[5];#chr:breakpoint fragment 2
    my $offset = $_[6]; #correct for overlap among fragments
    my $offset1 = $_[8];
    
    #calculate breakpoint coordinates using clipping order (NC=not clipped followed by clipped; CN=clipped followed by non-clipped) and offset
    #Four orientations possible...
    #1) NC <-> NC [2nd fragment is reverse complement]
    #2) CN <-> CN [2nd fragment is reverse complement]
    #3) NC <-> CN [both fragments forward orientation]
    #4) CN <-> NC [both fragments forward orientation]
    
    if ($clipOrder eq $clipOrder2){
        if($clipOrder eq "NC"){
            $bp_ += $offset1;
			$g1start = ($bp_ - $buffer);
			$g1end = $bp_;
			$bp2_ += $offset;
			$g2start = ($bp2_ - $buffer);
			$g2end = $bp2_;
        }else{
            $bp_ += $offset1;
			$g1start = $bp_ - 1;
            $g1end = $bp_ + ($buffer - 1);
            $bp2_ += $offset;
            $g2start = $bp2_ - 1;
            $g2end = $bp2_ + ($buffer - 1);
        }
    }else{
        if($clipOrder eq "NC"){
            $bp_ += $offset1;
            $g1start = ($bp_ - $buffer);
            $g1end = $bp_;
            $bp2_ += $offset;
            $g2end = $bp2_ + ($buffer - 1);
            $g2start = $bp2_ - 1;
        }else{
            $bp_ += $offset1;
            $g1start = $bp_ - 1;
            $g1end = $bp_ + ($buffer - 1);
            $bp2_ += $offset;
            $g2end = $bp2_;
            $g2start = ($bp2_ - $buffer);
        }
    }
    my $index1 = $bp;
    $index1 =~ s/:.*//g;
    $index1 = $index1 . ":" . $g1start . "-" . $g1end;
    my $index2 = $bp2;
    $index2 =~ s/:.*//g;
    $index2 = $index2 . ":" .$g2start . "-" . $g2end;
    if($g1start < 0 || $g2start < 0) { return ">$index1,$index2\nNNNNNN";}
    #extract genomic sequence for putative translocation and write to temporary file
    my $output = basename($bam);
    $output =~ s/bam/factera.tmp.fa/g;
    $output = "$OUTPUTDIR/$output";
    system("twoBitToFa -noMask $twobit:$index1,$index2 $output");
    
    #open temporary fasta file containing genomic sequences
    open(FILE, $output) or die $!;
    
    my $seqCount = 0;
    my $seq1 = "";
    my $seq2 = "";
    while(<FILE>){
        my $line = $_;
        chomp($line);
        if($line =~ m/>/){
			$seqCount++;
			next;
        }
        if($seqCount == 1){
			$seq1 = $seq1 . $line;
        }else {$seq2 = $seq2 . $line;}
    }
    
    close(FILE);
    #remove temporary file
    system("rm $output");
    
    #merge both sequences to create putative fusion gene sequence
    my $fusionseq = $seq1; 
    if($clipOrder eq $clipOrder2) { #need to flip fragment 2 into reverse complemen t(seq2)
        $seq2 = doRevComp($seq2);
    }
    if($clipOrder eq "NC") {$fusionseq = $fusionseq . $seq2;} #order=fragment 1, fragment 2
    else {$fusionseq = $seq2 . $fusionseq;} #order=fragment 2, fragment 1
    
    return ">$index1,$index2\n$fusionseq";
    
}

#===========================================================================================================       
#compare clipped read to non-clipped read using k-mer hashtable

sub compKmers{
    my $k = $_[3]; #size of k-mer
    my $sumk = 0; #count forward hits
    my $sumkrc = 0; #count reverse complement hits
    my $read1len = length($_[5]); #read 1 length ($_[5] is read1 seq)
    my $doOffset = $_[6];
    my $clip1 = $_[7]; #clip order, read 1
    my $clip2 = $_[8]; #clip order, read 2
    my $offset = 0;
    my $offsetrc = 0;
    my $firstmatch = 0;
    my $firstmatchrc = 0;

    for(my $itor = 0; $itor < length($_[2]) - $k; $itor++){
	my $km = substr($_[2], $itor, $k);

        if(exists($_[0]{$km})) {
            $sumk++;
            if($firstmatch == 0) {
                if($clip1 eq "CN" && $clip2 eq "NC"){
                    $offset = $itor - $_[0]{$km};
                }elsif($clip1 eq "NC" && $clip2 eq "CN") {
                    $offset = ($read1len - $_[0]{$km}) - (length($_[2]) - $itor);
                }
                $firstmatch = 1;
            }
        }
        if(exists($_[1]{$km})) {
            $sumkrc++;
            if($firstmatchrc == 0) {
                if($clip1 eq "NC" && $clip2 eq "NC"){
                    $offsetrc = $itor - $_[1]{$km};
                }elsif($clip1 eq "CN" && $clip2 eq "CN") {
                    $offsetrc = ($read1len - $_[1]{$km}) - (length($_[2]) - $itor);
		    #print "$offsetrc\t$read1len\t$itor\t" . $_[1]{$km} . "\t" . length($_[2]) . "\t$km\n";
                }
                $firstmatchrc = 1;
            }
        }
    }

    if($sumk >= max($_[4],($clipsize-$k))){
        my @tmp = ($offset,'F');
        return @tmp;
    }elsif($sumkrc >= max($_[4],($clipsize-$k))) {
        my @tmp = ($offsetrc,'RC');
        return @tmp;
    }
    else{return 9999999;}
}

#===========================================================================================================       
#convert read into kmers; store as hashtable

sub getKmers{
    my $k = $_[1]; #size of k
    my %kmers = (); #store k-mers

    for(my $itor = 0; $itor < length($_[0]) - $k; $itor++){
	my $km = substr($_[0], $itor, $k);
	if(!exists($kmers{$km})) {$kmers{$km} = $itor;}
    }
    return %kmers;
}                                                                                            

#=========================================================================================================== 
#parse cigar string                                                                                                           
sub parsecigar{
    my $cigar = $_[0];
    my $cigarItor = 0;
    my $var = substr($cigar,$cigarItor,1);
    my $num = "";
    #extract number and cigar code                                                                                            
    while($var =~ /^[+-]?\d+$/ && $cigarItor < length($cigar)){
        $cigarItor++;
        $num = $num . $var;
        $var = substr($cigar,$cigarItor,1);
    }

    $cigar = substr($cigar,($cigarItor+1));

    my @data;
    $data[0] = $num;
    $data[1] = $var;
    $data[2] = $cigar;

    return @data;
}


#===========================================================================================================
#return reverse complement of read
sub doRevComp
{
    my $rc = "";
    for(my $itor = length($_[0])-1; $itor >= 0; $itor--){
	my $c = substr($_[0], $itor, 1);
	if($c eq "A") {$c = "T";}
	elsif($c eq "T") {$c = "A";}
	elsif($c eq "C") {$c = "G";}
	elsif($c eq "G") {$c = "C";}
	$rc = $rc . $c;
    }
    return $rc;
}


#============================================================================================================
#given start coordinate of a read, find and return the closest gene------------------------------------------
sub findGene
{
    my $pos = $_[0];
    my $chr = $_[5];

    my $start_coor = BinSearch($pos, \&cmpFunc, \@{$_[1]});

    if((int $start_coor) - $start_coor != 0){
	my $p1 = $start_coor - 0.5;
	my $p2 = $start_coor + 0.5;
	my $pos1 = @{$_[1]}[$p1]; #-1
	my $pos2 = @{$_[1]}[$p2]; #-1

	if(abs($pos1 - $pos) < abs($pos2 - $pos)){
	    $start_coor = $pos1;
	}else {$start_coor = $pos2;}
    }else{
	$start_coor = @{$_[1]}[$start_coor];
    }
    my $start_coor_end = $getothercoor{"$chr $start_coor"};

    $_[2] = $start_coor;
    $_[3] = $start_coor_end;
    $_[4] = $coors2gene{"$chr $start_coor $start_coor_end"};

}


sub cmpFunc
{
    my ($index, $arrayRef, $target) = @_;
    my $item = $$arrayRef[$index];

    return $item <=> $target;
}

#binary search routine
sub BinSearch
{
    my ($target, $cmp) = @_;

    my $posmin = 0;
    my $posmax = $#{$_[2]};

    return -0.5 if &$cmp (0, \@{$_[2]}, $target) > 0;
    return $#{$_[2]} + 0.5 if &$cmp ($#{$_[2]}, \@{$_[2]}, $target) < 0;

while (1)
{
    my $mid = int (($posmin + $posmax) / 2);
    my $result = &$cmp ($mid, \@{$_[2]}, $target);
  
  if ($result < 0)
  {
      $posmin = $posmax, next if $mid == $posmin && $posmax != $posmin;
      return $mid + 0.5 if $mid == $posmin;
      $posmin = $mid;
  }
  elsif ($result > 0)
  {
      $posmax = $posmin, next if $mid == $posmax && $posmax != $posmin;
      return $mid - 0.5 if $mid == $posmax;
      $posmax = $mid;
  }
  else
  {
      return $mid;
  }
}
}

sub is_integer {
    local $_ = shift;
    # checks from perlfaq4
    return $] < 5.009002 unless defined;
    return 1 if (/^[+-]?\d+$/); # is a +/- integer
    0;
}

#get consensus sequence for soft-clipped reads for a given breakpoint
sub getConsensus
{
    my @seqarray = @{$_[0]{$_[1]}->{$_[2]}};
    my $sq = "";
    for(my $j = 0; $j < @seqarray; $j++) {

        my $mx_val = 0;
        my $mx_itor = 0;
        my $mx_char = "A";
        for(my $r = 0; $r < 4; $r++){
            if($seqarray[$j][$r] > $mx_val){
                $mx_val = $seqarray[$j][$r];
                $mx_itor = $r;
            }
        }
        next if $mx_val == 0;
        if($mx_itor == 1) {$mx_char = "T";}
        if($mx_itor == 2) {$mx_char = "C";}
        if($mx_itor == 3) {$mx_char = "G";}
    
        $sq = $sq . $mx_char;
    }
    return $sq;
}

#store consensus sequence for soft-clipped reads for a given breakpoint
sub getSeqArray
{
    my $i = $_[3];
    my $len = $_[4];
    my $j = $_[5];
    my $len2 = $_[6];
    
    my @seqarray = ();

    if(exists($_[7]{$_[1]}->{$_[2]})){
        @seqarray = @{$_[0]{$_[1]}->{$_[2]}};
    }
 
    for(; $i < $len; $i++){
        next if $i < 0;
        
        my $ch = substr($_[8], $j, 1);
        if($ch eq "A") {$seqarray[$i][0]++;}
        if($ch eq "T") {$seqarray[$i][1]++;}
        if($ch eq "C") {$seqarray[$i][2]++;}
        if($ch eq "G") {$seqarray[$i][3]++;}
        $j++;
        if($j >= $len2) {last;}
    }
    return @seqarray;

}


#perform bp correction by checking subsequence against reference
sub bp_correction
{
    my $read2 = $_[0];
    my @same = @{$_[1]};
    my $offsetbp2 = $_[2];
    my $bp2_ = $_[3];
    my $chr2 = $_[4];
    my $offsetbp1 = 0;
    
    my $index = "";
    my $r2seq = "";
    my @r2 = split(" ",$read2);
    if($same[0] > 0){
        #grab from read 2 downstream of cutpoint
        $r2seq = substr($r2[1], 0, $same[0]);
        $index = $chr2 . ":" . ($bp2_-1) . "-" . ($bp2_+$same[0]-1);
        if($bp2_-1 < 0) {return @data;}
    }else{
        #grab from read 2 upstream of cutpoint
        $r2seq = substr($r2[0], length($r2[0])+$same[0],-$same[0]);
        $index = $chr2 . ":" . ($bp2_+$same[0]) . "-" . $bp2_;
        if($bp2_+$same[0] < 0) {return @data;}
    }
    
    if(length($r2seq) < 3) {return ($offsetbp1, $offsetbp2);}
    
    #extract genomic sequence for putative translocation and write to temporary file
    my $output = basename($bam);
    $output =~ s/bam/factera.tmp.fa/g;
    $output = "$OUTPUTDIR/$output";
    system("twoBitToFa -noMask $twobit:$index $output");
    
    #open temporary fasta file containing genomic sequences
    open(FILE, $output) or die $!;
    
    my $seq1 = "";
    while(<FILE>){
        my $line = $_;
        chomp($line);
        next if($line =~ m/>/);
        $seq1 = $seq1 . $line;
    }
    
    close(FILE);
    #remove temporary file
    system("rm $output");

    if($r2seq eq $seq1){
        #adjust breakpoint 1
        if($same[1] eq 'F') {$offsetbp1 = -$offsetbp2;}
        else {$offsetbp1 = $offsetbp2;}
        $offsetbp2 = 0;
    }
    my @data = ($offsetbp1, $offsetbp2);
    return @data;
}



#check for 'chr' prefix on chromosome input argument; add if not present
sub addprefix
{
    if($disablechrprefix == 1) {return $_[0];}
    if(!($_[0] =~ m/chr/)){
        $_[0] = "chr" . $_[0];
    }else {$addchrprefix = 0;}
    
    return $_[0];
}


