#!/usr/bin/env perl
use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);
use Statistics::Basic qw(:all);

# Global variable
my ($help, $outDir, $line, $unit);

# Get Parameter
GetOptions(
	"h|help" => \$help,	
	"outDir=s" => \$outDir,
	"line=i" => \$line,
	"unit=i" => \$unit,
);

my $tmpDir = `pwd`;
chomp $tmpDir;
$outDir ||= "$tmpDir";
$line ||= 10000;
$unit ||= 50;

# Guide for program
my $guide_separator = "#" x 80;
my $program = basename(abs_path($0));
$program =~ s/\.pl//;
my $guide=<<INFO;
VER
	AUTHOR: xuxiangyang(xy_xu\@foxmail.com)
	NAME: $program
	PATH: $0
        VERSION: v0.1	2017-10-19
NOTE

USAGE
	$program <options> *.bam 
	$guide_separator Basic $guide_separator
	--help			print help information
	--outDir <s>		script out Dir, default "$outDir"
	--line <i>		how many number bam's lines used to count insert size, default $line
	--unit <i>		unit for count reads number, default $unit 

INFO

die $guide if (@ARGV == 0 || defined $help);

# Main
my @numbers;
my $count;

open BAM, "samtools view $ARGV[0]|" or die $!;
while (<BAM>) {
	$count++;
	last if($count > $line);
	my @arr = split;
	push @numbers, $arr[8] if($arr[8] > 0);
}
close BAM;

my $len = @numbers;
my $effective_rate = sprintf("%.2f%%", $len/$line*100);
my $median = median(@numbers);
my $mean = mean(@numbers);
my $stddev = stddev(@numbers);

my $sample = $ARGV[0];
$sample =~ s/.bam//i;
print "sample: $sample\n";
print "stat line: $line\n";
print "effective rate: $effective_rate\n";
print "median: $median\tmean: $mean\tstddev: $stddev\n";

my %region;
foreach my $key (@numbers) {
	my $value = int($key/$unit)*$unit;
	$region{$value}++;
}

my $file = basename $ARGV[0];
$file =~ s/.bam//i;

open RG, ">$outDir/$file.insertsize.xls" or die $!;
foreach my $range (sort {$a<=>$b} keys %region) {
	print RG "$range\t$region{$range}\n" if($range < $median*3);
}
close RG;

my $R=<<SCRIPT;

InsertData <- data.frame(read.table("$outDir/$file.insertsize.xls"))
png("$outDir/$file.insertsize.png", width = 1200, height = 500, res = 100)
barplot(InsertData[,2], width = 1, names.arg = InsertData[,1], xlab = "Insert Size", ylab = "Reads Pair Number", main = "Insert Size Evaluate", border = "blue", col = "blue", sub = "Stat for $line read pairs")
dev.off()

SCRIPT

open RS, ">$outDir/$file.insertsize.R" or die $!;
print RS $R;
close RS;

system("Rscript $outDir/$file.insertsize.R") == 0 or die $!;
