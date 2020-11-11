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
my ($help, $outDir, $segment, $region);
$outDir = $ENV{'PWD'};

# Get Parameter
GetOptions(
	"h|help" => \$help, 
	"s|segment" => \$segment,
	"outDir=s" => \$outDir,
	"region=s" => \$region,
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
	--segment               use segmentation calls (.cns) for copy number count, the program will automatically detect in the cnr directory
	--outDir <s>            script out Dir, default "$outDir"
	--region <s>            define which region will be draw, IGV-like format “chr:start-end”

INFO

die $guide if (@ARGV == 0 || defined $help);

# Main
system("echo \"$indent$0 Start at:`date '+%Y/%m/%d  %H:%M:%S'`$indent\"") == 0 || die $!;

foreach my $cnr_file (@ARGV) {
	my $prefix = basename($cnr_file);
	$prefix =~ s/.cnr//;
	my %cnr = ();
	my %cns = ();
	my @cn_x = ();
	my ($cnr_head, $cns_head, $cnr_index, $cns_index);
	open IN, $cnr_file or die $!;
	while (<IN>) {
		chomp;
		my @arr = split;
		if (/^chromosome/) {
			$cnr_head = "$_\tcn\tzscore\n";
			for (my $i = 0; $i < @arr; $i++) {
				if ($arr[$i] eq "log2") {
					$cnr_index = $i;
					last;
				}
			}
			next;
		}
		my ($cn, $zscore);
		my $log2 = $arr[$cnr_index];
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
		# set edge value for log2ratio and z-score
		$arr[$cnr_index] = -2 if ($arr[$cnr_index] < -2);
		$arr[$cnr_index] = 2 if ($arr[$cnr_index] > 2);
		$zscore = -3 if ($zscore < -3);
		$zscore = 3 if ($zscore > 3);
		if ($arr[0] =~ m/y/i) {
			#next if ($log2 < -2 or $log2 > 2);
			#next if ($zscore < -3 or $zscore > 3);
			$cn = $cn/2;
		}
		my $cnr_line = join "\t", @arr;
		$cnr{$arr[0]} .= "$cnr_line\t$cn\t$zscore\n";
		push @cn_x, $cn if ($arr[0] =~ m/x/i);
	}
	close IN;
	if (defined $segment) {
		my $cns_file = $cnr_file;
		$cns_file =~ s/.cnr$/.cns/;
		open SG, $cns_file or die $!;
		while (<SG>) {
			chomp;
			my @arr = split;
			if (/^chromosome/) {
				$cns_head = "$_\tcn\n";
				for (my $i = 0; $i < @arr; $i++) {
					if ($arr[$i] eq "log2") {
						$cns_index = $i;
						last;
					}
				}
				next;
			}
			my $cn;
			my $log2 = $arr[$cns_index];
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
			if ($arr[0] =~ m/y/i) {
				$cn = $cn/2;
			}
			$cns{$arr[0]} .= "$_\t$cn\n";
			$arr[1] = $arr[2];
			my $line = join "\t", @arr;
			$cns{$arr[0]} .= "$line\t$cn\n";
		}
		close SG;
	}
	my $sex = infer_sex(@cn_x); 
	#print "$cnr_file\t$sex\n";
	my @threads = ();
	foreach my $chr (sort {$a cmp $b} keys %cnr) {
		system("mkdir -p $outDir/$prefix.cnv.view") == 0 || die $!;
		open CNR, ">$outDir/$prefix.cnv.view/$prefix.$chr.cnr" or die $!;
		print CNR $cnr_head;
		print CNR $cnr{$chr}; 
		close CNR;
		if (defined $segment) {
		open CNS, ">$outDir/$prefix.cnv.view/$prefix.$chr.cns" or die $!;
		print CNS $cns_head;
		print CNS $cns{$chr}; 
		close CNS;
		}
		if ($region) {
			$region =~ m/^([^:]*)/;
			my $chr_region = $1; 
			push @threads, threads->create(\&draw_chr, "$outDir/$prefix.cnv.view/$prefix.$chr.cnr", $sex, $chr, "$outDir/$prefix.cnv.view", $segment, $region) if ($chr eq $chr_region); 
		} else {
			push @threads, threads->create(\&draw_chr, "$outDir/$prefix.cnv.view/$prefix.$chr.cnr", $sex, $chr, "$outDir/$prefix.cnv.view", $segment, $region); 
		}
	}			
	foreach (@threads) {
		$_->join();
	}
}

system("rm -rf $outDir/*.cnv.view/*.cnr $outDir/*.cnv.view/*.cns $outDir/*.cnv.view/*.R") == 0 || die $!;

system("echo \"$indent$0 End at:`date '+%Y/%m/%d  %H:%M:%S'`$indent\"") == 0 || die $!;

sub infer_sex {
	my $ave = average(@_);
	my $sex = "";
	if (1.6 <= $ave and $ave <= 2.4) {
		$sex = "Female";
	} elsif (0.6 <= $ave and $ave <= 1.4) {
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
	my $cnr_file = shift;
	my $sex = shift;
	my $chr = shift;
	my $dir = shift;
	my $seg = shift;
	my $reg = shift;
	my $pre = basename($cnr_file);
	$pre =~ s/.cnr//;
	$pre .= ".$reg" if (defined $reg);
	my $cns_file = $cnr_file;
	$cns_file =~ s/cnr$/cns/;
	my ($data, $cn_code);
	if (defined $seg) {
		$data=<<DATA;
cnr <- read.table("$cnr_file", header = T, sep = "\t", check.names = F)
cns <- read.table("$cns_file", header = T, sep = "\t", check.names = F)
pos <- cnr\$start
cns_pos <- cns\$start
cn <- cnr\$cn
cns_cn <- cns\$cn
DATA
		$cn_code=<<CNDE;	
kpLines(kp, chr = "$chr", x = cns_pos, y = cns_cn, r0 = 0.7, r1 = 1, ymin = 0, ymax = 4, col = "blue", cex = 0.5)
CNDE
	} else {
		$data=<<DATA;
cnr <- read.table("$cnr_file", header = T, sep = "\t", check.names = F)
pos <- cnr\$start
cn <- cnr\$cn
DATA
		$cn_code=<<CNDE;	
kpPoints(kp, chr = "$chr", x = pos, y = cn, r0 = 0.7, r1 = 1, ymin = 0, ymax = 4, col = "blue", cex = 0.4)
CNDE
	}
	my ($plotkt, $bnarg);
	if (defined $reg) {
		$plotkt = "zoom = \"$reg\"";
		my @tmp = split ":", $reg;
		my @tmp2 = split "-", $tmp[1];
		my $reg_len = $tmp2[1] - $tmp2[0];
		if ($reg_len < 100000) {
			$bnarg = "tick.dist = 10000, minor.tick.dist = 1000";
		} elsif (100000 <= $reg_len and $reg_len < 10000000) {
			$bnarg = "tick.dist = 1000000, minor.tick.dist = 100000";
		} else {
			$bnarg = "tick.dist = 10000000, minor.tick.dist = 1000000";
		}
	} else {
		$plotkt = "chromosomes = \"$chr\"";
		$bnarg = "tick.dist = 10000000, minor.tick.dist = 1000000";
	}
	my $Rscript=<<RS;
# load data
$data
zscore <- cnr\$zscore
log2ratio <- cnr\$log2

# set env
suppressMessages(library(karyoploteR))

pp <- getDefaultPlotParams(plot.type = 1)
pp\$bottommargin <- 50
pp\$topmargin <- 50
pp\$ideogramheight <- 20
pp\$data1height <- 400

# draw
png(file="$dir/$pre.png", width = 1200, height = 800, units = "px")

kp <- plotKaryotype(genome = "hg19", $plotkt, plot.type = 1, plot.params = pp, cex = 2)
kpAddBaseNumbers(kp, $bnarg, tick.len = 10, tick.col="red", minor.tick.len = 5, minor.tick.col = "gray", add.units = T, cex = 1)
kpAddCytobandLabels(kp, cex=1.2)
kpDataBackground(kp, r0 = 0, r1 = 1, color = 'white')

kpAxis(kp, r0 = 0.7, r1 = 1, ymin = 0, ymax = 4, cex = 1)
kpAddLabels(kp, labels = "Copy Number", r0 = 0.7, r1 = 1, cex = 1, label.margin = 0.035)
$cn_code

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

while (!is.null(dev.list()))  dev.off()
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

while (!is.null(dev.list()))  dev.off()
TL
	}
	open RO, ">$dir/$pre.R" or die $!;
	print RO $Rscript;
	close RO;
	system("Rscript $dir/$pre.R") == 0 || die $!;
	return 1;
}
