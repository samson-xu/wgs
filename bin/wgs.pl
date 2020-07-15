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

# File or tool path check
my $config = path_check("$Bin/config.txt");

# Global variable
my ($help, $stat, $target_region, $fastqc_help, $fastp_help, $backtrack_help, $mem_help, $mem2_help, %wgs_shell);
my $project = strftime("%Y%m%d",localtime());
my $interval_padding = 200;
my $workDir = $ENV{'PWD'};
my $ref = $config->{'hg19'};
my $step = '123456789';
my $thread = '35';
my $fastqc_arg = '';
my $fastp_arg = "--adapter_sequence AGATCGGAAGAG --adapter_sequence_r2 AGATCGGAAGAG -q 15 -u 40 -n 5 -l 50 -w $thread -d 3";
my $clean_fastq_split = 3;
my $align_way = 'mem2';
my $align_arg = '';

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
$indent 1. Quality control and check(FastQC) of input data(FASTQ).
$indent 2. Adapter cut and low quality sequence filter of fastq(Fastp).
$indent 3. Fastq Alignment and quality control of sequence alignment results.
$indent 4. SNP/InDel detection.
$indent 5. CNV detection.
$indent 6. Fusion gene detection.
$indent 7. SV detection.
$indent 8. Mitochondrial gene mutation detection.
$indent 9. Statistics of variation detection results

PARAMETER
$indent $0 [options] sample.lst

$parameter_separator Basic $parameter_separator 
$indent -help                        Print this guide information 
$indent -project <str>               Project name, default "$project"
$indent -target_region <str>         Target region bed files on the genome
$indent -interval_padding <i>        Amount of padding (in bp) to add to each interval, default "$interval_padding"
$indent -workDir <str>               Work directory, default "$workDir"
$indent -ref <str>                   Reference genome absolute path, default "$ref"
$indent -step <str>                  Set step for run, default "$step"
$indent -thread <i>                  Set the number of threads for the program to run, default "$thread"
$indent -stat                        Wether stat sample information, default not stat
$parameter_separator Filter $parameter_separator 
$indent -fastqc_help                 Print fastqc help information
$indent -fastqc_arg                  Fastqc argument setting, default "$fastqc_arg"
$indent -fastp_help                  Print fastp help information
$indent -fastp_arg                   Fastp argument setting, default "$fastp_arg"
$indent -clean_fastq_split           Specifies how many parts the output fastq is divided into, default "$clean_fastq_split"
$parameter_separator Align $parameter_separator 
$indent --align_way                  Select align algorithm, 'backtrack', 'mem' or 'mem2', default "$align_way"
$indent --backtrack_help             Print BWA-backtrack help information
$indent --mem_help                   Print BWA-mem help information
$indent --mem2_help                  Print BWA-mem2 help information
$indent --align_arg                  Align argument setting, this has to correspond to the align_way, default "$align_arg"

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
);

die $guide if ((@ARGV == 0 || defined $help) && !defined $fastqc_help && !defined $fastp_help && !defined $backtrack_help && !defined $mem_help && !defined $mem2_help);
die `$config->{fastqc} -h` if (defined $fastqc_help);
die `$config->{fastp}` if (defined $fastp_help);
die `$config->{bwa} aln` if (defined $backtrack_help);
die `$config->{bwa} mem` if (defined $mem_help);
die `$config->{mem2} mem` if (defined $mem2_help);

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
		if ($step =~ /1/) {
			my $fastqcDir = "$projectDir/$sampleId/00.fastqc";
			my $fastqc_shell = "$config->{fastqc} -t $thread -o $fastqcDir $fastqc_arg $fastq->[0] $fastq->[1]\n";
			$fastqc_shell .= "rm $fastqcDir/*.zip";
			write_shell($fastqc_shell, "$fastqcDir/$sampleId.fastqc.sh");
			$step1_shell .= "sh $fastqcDir/$sampleId.fastqc.sh >$fastqcDir/$sampleId.fastqc.sh.o 2>$fastqcDir/$sampleId.fastqc.sh.e";
		}
		# Fastq filter
		if ($step =~ /2/) {
			my $filterDir = "$projectDir/$sampleId/01.filter";
			my $filter_shell = "$config->{fastp} -i $fastq->[0] -o $filterDir/$sampleId.clean.1.fq.gz -I $fastq->[1] -O $filterDir/$sampleId.clean.2.fq.gz $fastp_arg -j $filterDir/$sampleId.fastq.json -h $filterDir/$sampleId.fastq.html -R '$sampleId fastq report'";
			write_shell($filter_shell, "$filterDir/$sampleId.filter.sh");
			$wgs_shell{$sampleId} .= $step1_shell . " &\n"if ($step =~ /1/);
			$wgs_shell{$sampleId} .= "sh $filterDir/$sampleId.filter.sh >$filterDir/$sampleId.filter.sh.o 2>$filterDir/$sampleId.filter.sh.e &\n";
			$wgs_shell{$sampleId} .= "\nwait\n\n";
			@{$sampleInfo{$sampleId}{'clean'}} = ("$filterDir/$sampleId.clean.1.fq.gz", "$filterDir/$sampleId.clean.2.fq.gz", $clean_fastq_split);
		}
		# Alignment
		if ($step =~ /3/) {
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

print STDERR "Because your input is alignment result, the program will automatically skip the first few steps!\n" if (!$fastq_label and $step =~ /1|2|3/);	
foreach my $sampleId (sort {$a cmp $b} keys %sampleInfo) {
	$sampleInfo{$sampleId}{'align'} = "$projectDir/$sampleId/02.align/$sampleId.final.bam" unless ($sampleInfo{$sampleId}{'align'});
	# SNP/InDel detection
	if ($step =~ /4/) {
		my $callDir = "$projectDir/$sampleId/03.variant";
		call_variant($config->{'gatk'}, $sampleInfo{$sampleId}{'align'}, $ref, $config->{'dict'}, $config->{'dbsnp'}, $config->{'mills'}, $config->{'axiom'}, $config->{'hapmap'}, $config->{'omni'}, $config->{'thousand'}, $callDir, $target_region);
		$wgs_shell{$sampleId} .= "sh $callDir/$sampleId.variant.sh >$callDir/$sampleId.variant.sh.o 2>$callDir/$sampleId.variant.sh.e &\n";
		$sampleInfo{$sampleId}{'variant'} = "$callDir/$sampleId.filter.vcf.gz"; 
	}
	# CNV detection
	if ($step =~ /5/) {

	}
	# Fusion gene detection
	if ($step =~ /6/) {

	}
	# SV detection
	if ($step =~ /7/) {

	}
	# Mitochondrial gene mutation detection
	if ($step =~ /8/) {

	}
	# Statistics of variation detection results
	if ($step =~ /9/) {

	}
}	


foreach my $sampleId (sort {$a cmp $b} keys %sampleInfo) {
	$wgs_shell{$sampleId} .= "\nwait\n\n";
	write_shell($wgs_shell{$sampleId}, "$projectDir/$sampleId/$sampleId.sh");
}

stat_log($sample_total, $Bin) if (defined $stat);