#!/usr/bin/env perl
#########################################################################
# File Name: qc_format.pl
# Author: samson-xu
# mail: xy_xu@foxmail.com
# Created Time: Wed 14 Oct 2020 04:41:46 PM CST
#########################################################################

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Data::Dumper;

my %hash = ();
while (<>) {
	chomp;
	my @arr = split /\t/;
	my $key = shift @arr;
	my $value = join "\t", @arr;
	$hash{$key} = $value;
}

my @bases = split /\t/, $hash{'Clean_bases'};
my @new_bases;
foreach my $base (@bases) {
	$base = sprintf("%.2f", $base/1000000);
	push @new_bases, $base;
}
$hash{'Clean_bases'} = join "\t", @new_bases;

my $output=<<CONTENT;
Sample\t$hash{'Sample'}
Clean_bases(Mb)\t$hash{'Clean_bases'}
Target_bases(bp)\t$hash{'Target_bases'}
Target_coverage_rate\t$hash{'Target_coverage_rate'}
Average_effective_depth_on_target(X)\t$hash{'Average_effective_depth_on_target'}
Fraction_target_covered_at_least_10x\t$hash{'Fraction_target_covered_at_least_10x'}
Fraction_target_covered_at_least_20x\t$hash{'Fraction_target_covered_at_least_20x'}
Fraction_target_covered_at_least_30x\t$hash{'Fraction_target_covered_at_least_30x'}
Q20_rate\t$hash{'Q20_rate'}
Q30_rate\t$hash{'Q30_rate'}
CONTENT

print $output;
