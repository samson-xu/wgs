package ReadsAlign;

use File::Basename;
use Data::Dumper;
use ScatterChr;
use WriteShell;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(reads_align);

sub reads_align {
	my $program = shift;
	my $samtools = shift;
	my $gatk = shift; 
	my $thread = shift;
	my $fastq = shift;
	my $ref = shift;
	my $dict = shift;
	my $dbsnp = shift;
	my $dbindel = shift;
	my $way = shift;
	my $arg = shift;
	my $outDir = shift;
	my ($bwa, $merge, $mark, $bqsr, $apply, $shell, $tmp_shell);
	#$thread = int($thread/$fastq->[2]);
	my $prefix = basename($fastq->[0]);
	$prefix =~ s/.clean.1.fq.gz//;
	my @prefix_list;
	if ($fastq->[2] == 1) {
		push @prefix_list, "";
	} else {
		for(my $i=1; $i<=$fastq->[2]; $i++) {
			push @prefix_list, sprintf("%03d.", $i);
		}
	}
	my $markdup_input;
	$bwa = "# BWA for alignment\n";
	foreach my $pre (@prefix_list) {
		my $dir = dirname($fastq->[0]);
		my $fq1 = basename($fastq->[0]);
		$fq1 = "$dir/$pre$fq1";
		my $fq2 = basename($fastq->[1]);
		$fq2 = "$dir/$pre$fq2";
		my $out_pre = "$pre$prefix";
		$markdup_input .= "--INPUT $outDir/$out_pre.merge.bam \\\n";
		if ($way eq 'backtrack') {
			$arg = "-t $thread" if ($arg eq '');
			$bwa .= "$program aln $arg $ref $fq1 > $outDir/$out_pre.aln_sa1.sai &\n";
			$bwa .= "$program aln $arg $ref $fq2 > $outDir/$out_pre.aln_sa2.sai &\n";
			$bwa .= "wait\n";
			$bwa .= "$program sampe $ref $outDir/$out_pre.aln_sa1.sai $outDir/$out_pre.aln_sa2.sai $fq1 $fq2 | $samtools view -1 --threads $thread - > $outDir/$out_pre.bam\n\n";
		} else {
			$arg = "-K 100000000 -t $thread -Y" if ($arg eq '');
			$bwa .= "$program mem $arg $ref $fq1 $fq2 | $samtools view -1 --threads $thread - > $outDir/$out_pre.bam\n";
		}

		$merge=<<FastqToSam;
# Convert FASTQ to read unmapped BAM
$gatk FastqToSam \\
-F1 $fq1 \\
-F2 $fq2 \\
-O $outDir/$out_pre.unmapped.bam \\
-SM $prefix \\
-PL ILLUMINA \\
-RG $prefix

FastqToSam
		$merge.=<<MergeBam;
# Merge original input uBAM file with BWA-aligned BAM file
$gatk MergeBamAlignment \\
--VALIDATION_STRINGENCY SILENT \\
--EXPECTED_ORIENTATIONS FR \\
--ATTRIBUTES_TO_RETAIN X0 \\
--ALIGNED_BAM $outDir/$out_pre.bam \\
--UNMAPPED_BAM $outDir/$out_pre.unmapped.bam \\
--OUTPUT $outDir/$out_pre.merge.bam \\
--REFERENCE_SEQUENCE $ref \\
--SORT_ORDER "unsorted" \\
--IS_BISULFITE_SEQUENCE false \\
--ALIGNED_READS_ONLY false \\
--CLIP_ADAPTERS false \\
--MAX_RECORDS_IN_RAM 2000000 \\
--ADD_MATE_CIGAR true \\
--MAX_INSERTIONS_OR_DELETIONS -1 \\
--PRIMARY_ALIGNMENT_STRATEGY MostDistant \\
--UNMAPPED_READ_STRATEGY COPY_TO_TAG \\
--ALIGNER_PROPER_PAIR_FLAGS true \\
--UNMAP_CONTAMINANT_READS true

rm $outDir/$out_pre.bam $outDir/$out_pre.unmapped.bam
MergeBam
	write_shell($merge, "$outDir/$out_pre.merge.sh");
	
	$tmp_shell .= "sh $outDir/$out_pre.merge.sh >$outDir/$out_pre.merge.sh.o 2>$outDir/$out_pre.merge.sh.e &\n";
	} 

	write_shell($bwa, "$outDir/$prefix.bwa.sh");
	$shell .= "sh $outDir/$prefix.bwa.sh >$outDir/$prefix.bwa.sh.o 2>$outDir/$prefix.bwa.sh.e\n";
	$shell .= $tmp_shell;
	
	$mark=<<MarkDuplicates;
# Mark duplicate reads to avoid counting non-independent observations
$gatk MarkDuplicates \\
$markdup_input--OUTPUT $outDir/$prefix.markdup.bam \\
--METRICS_FILE $outDir/$prefix.metrics \\
--VALIDATION_STRINGENCY SILENT \\
--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \\
--ASSUME_SORT_ORDER "queryname" \\
--CREATE_MD5_FILE false

rm $outDir/*.merge.bam 

MarkDuplicates

	$mark.=<<SortSam;
# Sort BAM file by coordinate order and fix tag values for NM and UQ
$gatk SortSam \\
--INPUT $outDir/$prefix.markdup.bam \\
--OUTPUT /dev/stdout \\
--SORT_ORDER "coordinate" \\
--CREATE_INDEX false \\
--CREATE_MD5_FILE false \\
| \\
$gatk SetNmMdAndUqTags \\
--INPUT /dev/stdin \\
--OUTPUT $outDir/$prefix.markdup.sort.bam \\
--CREATE_INDEX true \\
--CREATE_MD5_FILE true \\
--REFERENCE_SEQUENCE $ref

rm $outDir/$prefix.markdup.bam $outDir/$prefix.metrics

SortSam
	write_shell($mark, "$outDir/$prefix.mark.sh");		

	# Generate sets of intervals for scatter-gathering over chromosomes	
	my @intervals = @{scatter_chr($dict)};
	# add the unmapped sequences as a separate line to ensure that they are recalibrated as well
	my @intervals_with_unmapped = (@intervals, '-L unmapped');
	# print Dumper \@intervals_with_unmapped;
	my $groups;
	$bqsr = "#Generate Base Quality Score Recalibration (BQSR) model by interval\n";
	foreach my $interval (@intervals) {
		my @arr = split /\s+/, $interval;
		my $group = $arr[1];
		$groups .= "-I $outDir/$prefix.$group.recal_data.csv \\\n";
		$bqsr.=<<BaseRecalibrator; 
$gatk BaseRecalibrator \\
-R $ref \\
-I $outDir/$prefix.markdup.sort.bam \\
--use-original-qualities \\
-O $outDir/$prefix.$group.recal_data.csv \\
--known-sites $dbsnp \\
--known-sites $dbindel \\
$interval &

BaseRecalibrator

	}
	$bqsr .= "wait\n\n";
	$bqsr.=<<GatherBQSRReports;
# Combine multiple recalibration tables from scattered BaseRecalibrator runs
$gatk GatherBQSRReports \\
$groups-O $outDir/$prefix.recal_data.csv

GatherBQSRReports
	write_shell($bqsr, "$outDir/$prefix.bqsr.sh");		
	
	my $groups_with_unmapped;
	$apply = "#Apply Base Quality Score Recalibration (BQSR) model by interval\n";
	foreach my $interval (@intervals_with_unmapped) {
		my @arr = split /\s+/, $interval;
		my $group = $arr[1];
		$groups_with_unmapped .= "--INPUT $outDir/$prefix.$group.recal.bam \\\n";	
		$apply.=<<ApplyBQSR;
$gatk ApplyBQSR \\
-R $ref \\
-I $outDir/$prefix.markdup.sort.bam \\
-O $outDir/$prefix.$group.recal.bam \\
$interval \\
-bqsr $outDir/$prefix.recal_data.csv \\
--static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \\
--add-output-sam-program-record \\
--use-original-qualities &

ApplyBQSR
	}

	$apply .= "wait\n\n";
	$apply .= "rm $outDir/$prefix.markdup.sort.*\n\n";
	$apply.=<<GatherBamFiles;
# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
$gatk GatherBamFiles \\
$groups_with_unmapped--OUTPUT $outDir/$prefix.final.bam \\
--CREATE_INDEX true \\
--CREATE_MD5_FILE true

rm $outDir/*.recal*

GatherBamFiles
	write_shell($apply, "$outDir/$prefix.apply.sh");	
	$shell.=<<SHELL;
wait
sh $outDir/$prefix.mark.sh >$outDir/$prefix.mark.sh.o 2>$outDir/$prefix.mark.sh.e 
sh $outDir/$prefix.bqsr.sh >$outDir/$prefix.bqsr.sh.o 2>$outDir/$prefix.bqsr.sh.e 
sh $outDir/$prefix.apply.sh >$outDir/$prefix.apply.sh.o 2>$outDir/$prefix.apply.sh.e
SHELL
	write_shell($shell, "$outDir/$prefix.align.main.sh");
}

1;
