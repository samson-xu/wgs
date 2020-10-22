#!/usr/bin/env perl
#########################################################################
# File Name: gatk_joint_call.pl
# Author: samson-xu
# mail: xy_xu@foxmail.com
# Created Time: Wed 21 Oct 2020 03:50:29 PM CST
#########################################################################

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use Cwd 'abs_path';
use FindBin qw($Bin);
use Thread;

# Global variable
my ($help, $groupId, $outDir, $region, $reference, $dbsnp, $gatk);
$outDir = $ENV{'PWD'};;
$region = "wgs";
$reference = "/data/bioit/biodata/xuxy/pipeline/wgs/db/hg19.fa";
$dbsnp ="/data/bioit/biodata/xuxy/pipeline/wgs/db/dbsnp_138.hg19.vcf.gz";
$gatk = "/data/bioit/biodata/xuxy/pipeline/wgs/soft/gatk-4.1.8.0/gatk";

# Get Parameter
GetOptions(
	"h|help" => \$help, 
	"groupId=s" => \$groupId,
	"outDir=s" => \$outDir,
	"region=s" => \$region,
	"reference=s" => \$reference,
	"dbsnp=s" => \$dbsnp,
	"gatk=s" => \$gatk,
);

# Guide for program
my $guide_separator = "#" x 75;
my $indent = "#" x 20;
my $program = basename(abs_path($0));
$program =~ s/\.pl//;
my $guide=<<INFO;
VER
	AUTHOR: xuxiangyang(xy_xu\@foxmail.com)
	NAME: $program
	PATH: $0
	VERSION: v1.0   2020-10-21
NOTE

USAGE
	$program <options> *.g.vcf.gz 
	$guide_separator Basic $guide_separator
	--help                  print help information
	--groupId               group id for gatk joint call, default automatic recognition
	--outDir <s>            script out Dir, default "$outDir"
	--region <s>            set the target interval(bed format) of genotyping, wgs dosen't need to be set
	--reference <s>         reference genome path, default "$reference"
	--dbsnp <s>             dbsnp path, default "$dbsnp"
	--gatk <s>              gatk path, default "$gatk"

INFO

die $guide if (@ARGV == 0 || defined $help);

# Main
system("echo \"$indent$0 Start at:`date '+%Y/%m/%d  %H:%M:%S'`$indent\"") == 0 || die $!;
my @files = map {"-V $_ \\"} @ARGV;
my $input = join "\n", @files;
# auto match family id
for (my $i=0; $i<@ARGV; $i++) {
	unless (-e "$ARGV[$i].tbi") {
		print STDERR "Please make index for $ARGV[$i]!\n";
		exit;
	}
	$ARGV[$i] =~ s/.gz$//;
	$ARGV[$i] =~ s/.g.vcf$//;
}
$groupId ||= list_com_str(@ARGV);
my @pool = @{scatter_normal_chr($region, $outDir)};
# parallel GenotypeGVCFs
my @threads; 
foreach (@pool) {
	push @threads, threads->create(\&GenotypeGVCFs, $input, $_, $region, $reference, $dbsnp, $outDir);
} 
foreach (@threads) {
	$_->join();
}
# GatherVcfs
my @gather_inputs = map {"-I $outDir/chr$_.vcf.gz \\"} (1..22, 'X', 'Y');
my $gather_input = join "\n", @gather_inputs;
my $gather_shell=<<GAS;
$gatk --java-options -Xms6g \\
GatherVcfs \\
$gather_input
-O $outDir/$groupId.vcf.gz
GAS
#print $gather_shell;
system("$gather_shell") == 0 || die $!;
system("echo \"$indent$0 End at:`date '+%Y/%m/%d  %H:%M:%S'`$indent\"") == 0 || die $!;
system("rm -rf $outDir/*chr* $outDir/bed_split") == 0 || die $!;


sub max_mutual_str {
	my ($min,$max) = @_;
	($min,$max) = ($max,$min) if (length($min) > length($max));
	my $len = length($min);
	foreach my $rest (0..$len-1) {
		foreach my $pos (0..$rest) {
			my $str = substr($min, $pos, ($len-$rest));
			return $str if $max =~ /$str/;
		}
	}
	return undef;
}

sub list_com_str {
	my @arr = @_;
	my $first = shift @arr;
	my $second = shift @arr;
	my $common_str = max_mutual_str($first, $second);
	if (@arr >0) {
		foreach my $str (@arr) {
			$common_str = max_mutual_str($common_str, $str);
		}
	}
	return $common_str;
}

sub scatter_normal_chr {
	my $file = shift;
	my $outDir = shift;
	$outDir = "$outDir/bed_split";
	my @intervals;
	my @chrs = map {"chr$_"} (1..22, 'X', 'Y');
	if ($file =~ /bed$/i) {
		my %chr_region;
		open F, $file or die $!;
		while (<F>) {
				next if (/#/);
				next if (/^\s*$/);
				chomp;
				my @arr = split /\s+/;
				$chr_region{$arr[0]} .= "$_\n";
		}
		close F;
		system("mkdir -p $outDir") == 0 || die $!;
		foreach my $chr (keys %chr_region) {
				open OUT, ">$outDir/$chr.bed" or die $!;
				print OUT $chr_region{$chr};
				close OUT;
		}
		foreach my $chr (@chrs) {
				push @intervals, "$outDir/$chr.bed" if (-e "$outDir/$chr.bed");
		}
	} else {
		@intervals = @chrs;	
	}
	return \@intervals;
}

sub GenotypeGVCFs {
	my $input = shift;
	my $chr = shift;
	my $target = shift;
	my $ref = shift;
	my $dbsnp = shift;
	my $outDir = shift;
	my $prefix = basename($chr); 
	$prefix =~ s/.bed//;
	my $interval_padding = "";
	$interval_padding = "-ip 100" if ($target =~ m/.bed$/i);
	my $genotype_shell=<<GTS;
$gatk --java-options -Xmx5g \\
GenomicsDBImport \\
$input
--genomicsdb-workspace-path $outDir/genomicsdb.$prefix \\
-L $chr $interval_padding \\
--reader-threads 5 \\
--overwrite-existing-genomicsdb-workspace \\
--merge-input-intervals \\
--consolidate

$gatk --java-options -Xmx5g \\
GenotypeGVCFs \\
-R $ref \\
-O $outDir/$prefix.vcf.gz \\
-D $dbsnp \\
-G AS_StandardAnnotation \\
--only-output-calls-starting-in-intervals \\
-V gendb://$outDir/genomicsdb.$prefix \\
-L $chr $interval_padding \\
--merge-input-intervals


GTS
	#print $genotype_shell;
	system("$genotype_shell") == 0 || die $!;
}
