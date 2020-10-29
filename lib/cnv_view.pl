#!/usr/bin/env perl
#########################################################################
# File Name: cnv_view.pl
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
my ($help, $outDir);
$outDir = $ENV{'PWD'};;

# Get Parameter
GetOptions(
	"h|help" => \$help, 
	"outDir=s" => \$outDir,
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
	VERSION: v1.0   2020-10-28
NOTE

USAGE
	$program <options> *.cnr 
	$guide_separator Basic $guide_separator
	--help                  print help information
	--outDir <s>            script out Dir, default "$outDir"

INFO

die $guide if (@ARGV == 0 || defined $help);

# Main
system("echo \"$indent$0 Start at:`date '+%Y/%m/%d  %H:%M:%S'`$indent\"") == 0 || die $!;
foreach my $file (@ARGV) {
	my $prefix = basename($file);
	$prefix =~ s/.cnr//;
	my %cnr = ();
	my @cn_x = ();
	my ($head, $index);
	open IN, $file or die $!;
	while (<IN>) {
		chomp;
		my @arr = split;
		if (/^chromosome/) {
			$head = "$_\tcn\tzscore\n";
			for (my $i = 0; $i < @arr; $i++) {
				if ($arr[$i] eq "log2") {
					$index = $i;
					last;
				}
			}
			next;
		}
		my ($cn, $zscore);
		my $log2 = $arr[$index];
		if ($log2 <= -1.1) {
			$cn = 0;
		} elsif ( -1.1 < $log2 and $log2 <= -0.4) {
			$cn = 1;	
		} elsif (-0.4 < $log2 and $log2 <= 0.3) {
			$cn = 2;
		} elsif (0.3 < $log2 and $log2 <= 0.7) {
			$cn = 3;
		} else {
			$cn = 4;
		}
		if ($prefix =~ m/ref.cnn/i) {
			$zscore = sprintf("%.4f", $log2/(1-$arr[-1]));
		} else {
			$zscore = sprintf("%.4f", $log2/sqrt(1-$arr[-1]));
		}
		if ($arr[0] =~ m/y/i) {
			next if ($log2 < -2 or $log2 > 2);
			next if ($zscore < -3 or $zscore > 3);
			$cn = $cn/2;
		}
		$cnr{$arr[0]} .= "$_\t$cn\t$zscore\n";
		push @cn_x, $cn if ($arr[0] =~ m/x/i);
	}
	close IN;
	my $sex = infer_sex(@cn_x); 
	#print "$file\t$sex\n";
	my @threads = ();
	foreach my $chr (sort {$a cmp $b} keys %cnr) {
		system("mkdir -p $outDir/$prefix.cnv.view") == 0 || die $!;
		open OUT, ">$outDir/$prefix.cnv.view/$prefix.$chr.cnr" or die $!;
		print OUT $head;
		print OUT $cnr{$chr}; 
		close OUT;
		push @threads, threads->create(\&draw_chr, "$outDir/$prefix.cnv.view/$prefix.$chr.cnr", $sex, $chr, "$outDir/$prefix.cnv.view"); 
	}			
	foreach (@threads) {
		$_->join();
	}
}
system("echo \"$indent$0 End at:`date '+%Y/%m/%d  %H:%M:%S'`$indent\"") == 0 || die $!;

sub infer_sex {
	my $ave = average(@_);
	my $sex = "";
	if (1.8 <= $ave and $ave <= 2.2) {
		$sex = "Female";
	} elsif (0.8 <= $ave and $ave <= 1.2) {
		$sex = "Male";
	} else {
		$sex = "Unknow";
	}
	return $sex;	
}

sub average {
	my $n = scalar(@_);
	my $sum = 0;
	foreach my $item (@_) {
		 $sum += $item;
	}
	my $ave = $sum / $n;
	return $ave;	
}

sub draw_chr {
	my $file = shift;
	my $sex = shift;
	my $chr = shift;
	my $dir = shift;
	my $pre = basename($file);
	$pre =~ s/.cnr//;
	my $Rscript=<<RS;
# load data
cnr <- read.table("$file", header = T, sep = "\t", check.names = F)
pos <- cnr\$start
cn <- cnr\$cn
zscore <- cnr\$zscore
log2ratio <- cnr\$log2

# set env
library(karyoploteR)

pp <- getDefaultPlotParams(plot.type = 1)
pp\$bottommargin <- 50
pp\$topmargin <- 50
pp\$ideogramheight <- 20
pp\$data1height <- 400

# draw
png(file="$dir/$pre.png", width = 1200, height = 800, units = "px")

kp <- plotKaryotype(genome = "hg19", chromosomes = "$chr", plot.type = 1, plot.params = pp, cex = 2)
kpAddBaseNumbers(kp, add.units = T, cex = 1)
kpDataBackground(kp, r0 = 0, r1 = 1, color = 'white')

kpAxis(kp, r0 = 0.7, r1 = 1, ymin = 0, ymax = 4, cex = 1)
kpAddLabels(kp, labels = "Copy Number", r0 = 0.7, r1 = 1, cex = 1, label.margin = 0.035)
kpPoints(kp, chr = "$chr", x = pos, y = cn, r0 = 0.7, r1 = 1, ymin = 0, ymax = 4, col = "blue", cex = 0.4)

RS
	if ($sex eq "Male" and $chr =~ m/x/i) {
	$Rscript.=<<TL;
kpAxis(kp, r0 = 0, r1 = 0.3, ymin = -5, ymax = 1, cex = 1)
kpAddLabels(kp, labels="Z Score", r0 = 0, r1 = 0.3, cex = 1, label.margin = 0.035)
kpPoints(kp, chr = "$chr", x = pos, y = zscore, r0 = 0, r1 = 0.3, ymin = -5, ymax = 1, col = "gray", cex = 0.5)

kpAxis(kp, r0 = 0.35, r1 = 0.65, ymin = -3, ymax = 1, cex = 1)
kpAddLabels(kp, labels = "Log2Ratio", r0 = 0.35, r1 = 0.65, cex = 1, label.margin = 0.035)
kpPoints(kp, chr = "$chr", x = pos, y = log2ratio, r0 = 0.35, r1 = 0.65, ymin = -3, ymax = 1, col = "green", cex = 0.5)
kpAbline(kp, h = -0.4, r0 = 0.35, r1 = 0.65, col = "gray", ymin = -3, ymax = 1, lwd = 1, lty = 2)
kpAbline(kp, h = -1.1, r0 = 0.35, r1 = 0.65, col = "gray", ymin = -3, ymax = 1, lwd = 1, lty = 2)

dev.off()
TL
	} else {
	$Rscript.=<<TL;
kpAxis(kp, r0 = 0, r1 = 0.3, ymin = -3, ymax = 3, cex = 1)
kpAddLabels(kp, labels="Z Score", r0 = 0, r1 = 0.3, cex = 1, label.margin = 0.035)
kpPoints(kp, chr = "$chr", x = pos, y = zscore, r0 = 0, r1 = 0.3, ymin = -3, ymax = 3, col = "gray", cex = 0.5)

kpAxis(kp, r0 = 0.35, r1 = 0.65, ymin = -2, ymax = 2, cex = 1)
kpAddLabels(kp, labels = "Log2Ratio", r0 = 0.35, r1 = 0.65, cex = 1, label.margin = 0.035)
kpPoints(kp, chr = "$chr", x = pos, y = log2ratio, r0 = 0.35, r1 = 0.65, ymin = -2, ymax = 2, col = "green", cex = 0.5)
kpAbline(kp, h = 0.3, r0 = 0.35, r1 = 0.65, col = "gray", ymin = -2, ymax = 2, lwd = 1, lty = 2)
kpAbline(kp, h = -0.4, r0 = 0.35, r1 = 0.65, col = "gray", ymin = -2, ymax = 2, lwd = 1, lty = 2)

dev.off()
TL
	}
	open RO, ">$dir/$pre.R" or die $!;
	print RO $Rscript;
	close RO;
	system("Rscript $dir/$pre.R") == 0 || die $!;
	system("rm $dir/$pre.R $dir/$pre.cnr") == 0 || die $!;
	return 1;
}
