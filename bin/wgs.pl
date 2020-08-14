#!/usr/bin/env perl
################################################################################
# File Name: wgs.pl
# Author: samson-xu
# mail: xy_xu@foxmail.com
# Created Time: Mon 22 Jun 2020 02:57:53 PM CST
################################################################################

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use POSIX;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use ConfigParse;
use SampleStat;
use WriteShell;
use ReadsAlign;
use VariantCall;
use SV;

# File or tool path check
my $config = path_check("$Bin/config.txt");

# Global variable
my ($help, $stat, $fastqc_help, $fastp_help, $backtrack_help, $mem_help, $mem2_help, $fusion_help, %wgs_shell, $main_shell);
my $project = strftime("%Y%m%d",localtime());
my $target_region = '';
my $interval_padding = 0;
my $workDir = $ENV{'PWD'};
my $ref = $config->{'hg19'};
my $step = 'cfbvnusmt';
my $thread = '35';
my $fastqc_arg = '';
my $fastp_arg = "--detect_adapter_for_pe -q 15 -u 40 -n 5 -l 50 -w $thread -d 3";
my $clean_fastq_split = 3;
my $align_way = 'mem2';
my $align_arg = '';
my $fusion_arg = '';

# Guide
my $guide_separator = "=" x 150;
my $indent = " " x 15;
my $parameter_separator = "*" x 70;
my $guide=<<INFO;
$guide_separator
=
=$indent$indent$indent NAME: Whole Genome Sequencing Analysis Pipeline(wgs) 
=$indent$indent$indent AUTHOR: xuxiangyang(xy_xu\@foxmail.com) 
=$indent$indent$indent VERSION: v0.1	2020/06/24                             
=                                          
$guide_separator


FUNCTIONS
$indent c. Quality control and check(FastQC) of input data(FASTQ).
$indent f. Adapter cut and low quality sequence filter of fastq(Fastp).
$indent b. Fastq Alignment and quality control of sequence alignment results.
$indent v. SNP/InDel detection.
$indent n. CNV detection.
$indent u. Fusion gene detection.
$indent s. SV detection.
$indent m. Mitochondrial gene mutation detection.
$indent t. Statistics of variation detection results

PARAMETER
$indent $0 [options] sample.lst

$parameter_separator Basic $parameter_separator 
$indent --help                        Print this guide information 
$indent --project <str>               Project name, default "$project"
$indent --target_region <str>         Target region bed files on the genome
$indent --interval_padding <i>        Amount of padding (in bp) to add to each interval, default "$interval_padding"
$indent --workDir <str>               Work directory, default "$workDir"
$indent --ref <str>                   Reference genome absolute path, default "$ref"
$indent --step <str>                  Set step for run, default "$step"
$indent --thread <i>                  Set the number of threads for the program to run, default "$thread"
$indent --stat                        Wether stat sample information, default not stat
$parameter_separator Filter $parameter_separator 
$indent --fastqc_help                 Print fastqc help information
$indent --fastqc_arg                  Fastqc argument setting, default "$fastqc_arg"
$indent --fastp_help                  Print fastp help information
$indent --fastp_arg                   Fastp argument setting, default "$fastp_arg"
$indent --clean_fastq_split           Specifies how many parts the output fastq is divided into, default "$clean_fastq_split"
$parameter_separator Align $parameter_separator 
$indent --align_way                  Select align algorithm, 'backtrack', 'mem' or 'mem2', default "$align_way"
$indent --backtrack_help             Print BWA-backtrack help information
$indent --mem_help                   Print BWA-mem help information
$indent --mem2_help                  Print BWA-mem2 help information
$indent --align_arg                  Align argument setting, this has to correspond to the align_way, default "$align_arg"
$parameter_separator Fusion $parameter_separator 
$indent --fusion_help                Print fusion help information 
$indent --fusion_arg                 Fusion argument setting, default "$fusion_arg"

NOTE
$indent 1. Fastq quality system should be phred 33
$indent 2. If input is fastq, sample.lst format like: SampleId    fq1    fq2, if input is bam, sample.lst format like: SampleId    bam
$indent 3. If you specify a target it will run as the WES, otherwise it will run the WGS
$guide_separator

INFO

# Parameter
GetOptions(
	"h|help" => \$help,
	"project=s" => \$project,
	"target_region=s" => \$target_region,
	"interval_padding=i" => \$interval_padding,
	"workDir=s" => \$workDir,
	"ref=s" => \$ref,
	"step=s" => \$step,
	"thread=i" => \$thread,
	"stat" => \$stat,
	"fastqc_help" => \$fastqc_help,
	"fastqc_arg=s" => \$fastqc_arg,
	"fastp_help" => \$fastp_help,
	"fastp_arg=s" => \$fastp_arg,
	"clean_fastq_split=i" => \$clean_fastq_split,
	"align_way=s" => \$align_way, 
	"backtrack_help" => \$backtrack_help,
	"mem_help" => \$mem_help,
	"mem2_help" => \$mem2_help,
	"fusion_help" => \$fusion_help,
	"fusion_arg=s" => \$fusion_arg,
);

#die $guide if ((@ARGV == 0 || defined $help) && !defined $fastqc_help && !defined $fastp_help && !defined $backtrack_help && !defined $mem_help && !defined $mem2_help && !defined $fusion_help);
if (@ARGV == 0) {
	die `$config->{fastqc} -h` . "\n" if (defined $fastqc_help);
	die `$config->{fastp}` . "\n" if (defined $fastp_help);
	die `$config->{bwa} aln` . "\n" if (defined $backtrack_help);
	die `$config->{bwa} mem` . "\n" if (defined $mem_help);
	die `$config->{mem2} mem` . "\n" if (defined $mem2_help);
	die `$config->{fusion}`."\n" if (defined $fusion_help);
}
die $guide if (@ARGV == 0 || defined $help);

$fastp_arg .= " --split $clean_fastq_split" if ($clean_fastq_split and $clean_fastq_split >=2);

# Main
my $projectDir = "$workDir/$project";
my $sample_file = shift;
system("mkdir -p $projectDir") == 0 || die $!;
system("cp $sample_file $projectDir") == 0 || die $!;
# load sample info
my %sampleInfo;
my $fastq_label = `grep 'gz\$' $sample_file`;
open SF, $sample_file or die $!;
while (<SF>) {
	next if (/#/);
	chomp;
	my @arr = split /\s+/;
	if ($fastq_label) {
		@{$sampleInfo{$arr[0]}{'fastq'}} = ($arr[1], $arr[2]);
	} else {
		$sampleInfo{$arr[0]}{'align'} = $arr[1];
	}
}
close SF;
#print Dumper \%sampleInfo;
my $sample_total = keys %sampleInfo;

my ($step1_shell);

if ($fastq_label) {
	foreach my $sampleId (sort {$a cmp $b} keys %sampleInfo) {
		# Fastq quality control
		my $fastq = $sampleInfo{$sampleId}{'fastq'};
		if ($step =~ /c/) {
			my $fastqcDir = "$projectDir/$sampleId/00.fastqc";
			my $fastqc_shell = "$config->{fastqc} -t $thread -o $fastqcDir $fastqc_arg $fastq->[0] $fastq->[1]\n";
			$fastqc_shell .= "rm $fastqcDir/*.zip";
			write_shell($fastqc_shell, "$fastqcDir/$sampleId.fastqc.sh");
			$step1_shell = "sh $fastqcDir/$sampleId.fastqc.sh >$fastqcDir/$sampleId.fastqc.sh.o 2>$fastqcDir/$sampleId.fastqc.sh.e";
		}
		# Fastq filter
		if ($step =~ /f/) {
			my $filterDir = "$projectDir/$sampleId/01.filter";
			my $filter_shell = "$config->{fastp} -i $fastq->[0] -o $filterDir/$sampleId.clean.1.fq.gz -I $fastq->[1] -O $filterDir/$sampleId.clean.2.fq.gz $fastp_arg -j $filterDir/$sampleId.fastq.json -h $filterDir/$sampleId.fastq.html -R '$sampleId fastq report'\n";
			$filter_shell .= "perl -I '$Bin/../lib' -MReadsStat -e \"reads_stat('$filterDir/$sampleId.fastq.json')\"\n";
			write_shell($filter_shell, "$filterDir/$sampleId.filter.sh");
			$wgs_shell{$sampleId} .= $step1_shell . " &\n"if ($step =~ /c/);
			$wgs_shell{$sampleId} .= "sh $filterDir/$sampleId.filter.sh >$filterDir/$sampleId.filter.sh.o 2>$filterDir/$sampleId.filter.sh.e &\n";
			$wgs_shell{$sampleId} .= "\nwait\n\n";
			@{$sampleInfo{$sampleId}{'clean'}} = ("$filterDir/$sampleId.clean.1.fq.gz", "$filterDir/$sampleId.clean.2.fq.gz", $clean_fastq_split);
		}
		# Alignment
		if ($step =~ /b/) {
			my $alignDir = "$projectDir/$sampleId/02.align";
			@{$sampleInfo{$sampleId}{'clean'}} = ($fastq->[0], $fastq->[1], 1) unless ($sampleInfo{$sampleId}{'clean'});
			my $align_program = $config->{'mem2'};
			$align_program = $config->{'bwa'} if ($align_way ne 'mem2');
			reads_align($align_program, $config->{'samtools'}, $config->{'gatk'}, $thread, $sampleInfo{$sampleId}{'clean'},
                        $ref, $config->{'dict'}, $config->{'dbsnp'}, $config->{'mills'}, $align_way, $align_arg, $alignDir);
			$wgs_shell{$sampleId} .= "sh $alignDir/$sampleId.align.sh >$alignDir/$sampleId.align.sh.o 2>$alignDir/$sampleId.align.sh.e\n";
			$sampleInfo{$sampleId}{'align'} = "$alignDir/$sampleId.final.bam"; 
		}	
	}
}

print STDERR "Because your input is alignment result, the program will automatically skip the first few steps!\n" if (!$fastq_label and $step =~ /c|f|b/);	
foreach my $sampleId (sort {$a cmp $b} keys %sampleInfo) {
	$sampleInfo{$sampleId}{'align'} = "$projectDir/$sampleId/02.align/$sampleId.final.bam" unless ($sampleInfo{$sampleId}{'align'});
	# SNP/InDel detection
	if ($step =~ /v/) {
		my $callDir = "$projectDir/$sampleId/03.variant";
		call_variant($config->{'gatk'}, $sampleInfo{$sampleId}{'align'}, $ref, $config->{'dict'}, $config->{'dbsnp'}, $config->{'mills'}, $config->{'axiom'}, $config->{'hapmap'}, $config->{'omni'}, $config->{'thousand'}, $callDir, $target_region, $interval_padding);
		$wgs_shell{$sampleId} .= "sh $callDir/$sampleId.variant.sh >$callDir/$sampleId.variant.sh.o 2>$callDir/$sampleId.variant.sh.e &\n";
		$wgs_shell{$sampleId} .= "perl -I '$Bin/../lib' -MBamQc -e \"bam_qc('$sampleInfo{$sampleId}{'align'}', '$config->{'samtools'}', '$target_region', '$config->{'bedtools'}', $thread)\" &\n";
		$wgs_shell{$sampleId} .= "\nwait\n\n";
		$sampleInfo{$sampleId}{'variant'} = "$callDir/$sampleId.filter.vcf.gz"; 
	}
	# CNV detection
	if ($step =~ /n/) {

	}
	# Fusion gene detection
	if ($step =~ /u/) {
		my $fusionDir = "$projectDir/$sampleId/04.sv/fusion";
		my $fusion_shell=<<FUSION;
# Fusion detect
$config->{'fusion'} \\
-o $fusionDir \\
-A $config->{'samtools'} \\
-T $config->{'twoBitToFa'} \\
-B $config->{'blastn'} \\
-M $config->{'makeblastdb'} \\
-p $thread \\
-F $fusion_arg \\
$sampleInfo{$sampleId}{'align'} \\
$config->{'interGene'} \\
$config->{'2bit'} $target_region
FUSION
		write_shell($fusion_shell, "$fusionDir/$sampleId.fusion.sh");
		$wgs_shell{$sampleId} .= "sh $fusionDir/$sampleId.fusion.sh >$fusionDir/$sampleId.fusion.sh.o 2>$fusionDir/$sampleId.fusion.sh.e\n";
		$sampleInfo{$sampleId}{'fusion'} = "$fusionDir/$sampleId.fusions.txt";
	}
	# SV detection
	if ($step =~ /s/) {
		my $mantaDir = "$projectDir/$sampleId/04.sv/manta";
		my $region = $config->{'dict'}; 
		$region = $target_region if ($target_region);
		manta($region, $config->{'manta'}, $config->{'convertInversion'}, $config->{'samtools'}, $config->{'AnnotSV'}, $sampleInfo{$sampleId}{'align'}, $config->{'hg19'}, $mantaDir, $thread, $config->{'bgzip'}, $config->{'tabix'});
		$wgs_shell{$sampleId} .= "sh $mantaDir/$sampleId.manta.sh >$mantaDir/$sampleId.manta.sh.o 2>$mantaDir/$sampleId.manta.sh.e\n";
		$sampleInfo{$sampleId}{'manta'} = "$mantaDir/results/variants/diploidSV.vcf.gz"; 
	}
	# Mitochondrial gene mutation detection
	if ($step =~ /m/) {

	}
	# Statistics of variation detection results
	if ($step =~ /t/) {

	}
}	

$main_shell = "# Run wgs pipeline for all samples\n";

foreach my $sampleId (sort {$a cmp $b} keys %sampleInfo) {
	write_shell($wgs_shell{$sampleId}, "$projectDir/$sampleId/$sampleId.sh");
	$main_shell .= "sh $projectDir/$sampleId/$sampleId.sh >$projectDir/$sampleId/$sampleId.sh.o 2>$projectDir/$sampleId/$sampleId.sh.e\n";
}

$main_shell.=<<Merge;
paste $projectDir/*/01.filter/*.fq.stat.txt | awk '{for(i=3; i<=NF; i+=2){\$i=""}; print \$0}' | sed "s/\\s\\+/\\t/g" > $projectDir/sample.fq.stat.xls
paste $projectDir/*/02.align/*.bam.stat.txt | awk '{for(i=3; i<=NF; i+=2){\$i=""}; print \$0}' | sed "s/\\s\\+/\\t/g" > $projectDir/sample.bam.stat.xls
cat $projectDir/sample.fq.stat.xls $projectDir/sample.bam.stat.xls | sed '7d' | sed '8,10d'> $projectDir/sample.stat.xls 
Merge

write_shell($main_shell, "$projectDir/main.sh");

stat_log($sample_total, $Bin) if (defined $stat);
