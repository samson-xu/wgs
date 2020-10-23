#!/usr/bin/env perl
#########################################################################
# File Name: gvcf_sex.pl
# Author: samson-xu
# mail: xy_xu@foxmail.com
# Created Time: Fri 23 Oct 2020 10:01:36 AM CST
#########################################################################

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use Thread;

if (@ARGV == 0) {
	print STDERR "$0 *.g.vcf\n";
	exit;
}

# parallel sex_infer
my @threads;
foreach (@ARGV) {
    push @threads, threads->create(\&sex_infer, $_);
}
foreach (@threads) {
    $_->join();
}

#sex_infer($ARGV[0]);

sub sex_infer {
	my $vcf = shift;
	my (%hash, %sample_index);
	if ($vcf =~ m/gz$/i) {
		open VCF, "gzip -dc $vcf |" or die $!;
	} else {
		open VCF, $vcf or die $!;
	}
	while (<VCF>) {
		next if (/^##/);
		chomp;
		my @arr = split /\t/;
		if (/^#CHROM/) {
			for (my $i=9; $i<@arr; $i++) {
				$sample_index{$i} = $arr[$i];
			}
		}
		my @alts = split /,/, $arr[4];
		next if (@alts > 2);
		next unless (/^chrX/i);
		next unless (/GT:AD:DP/ or /GT:DP/);
		foreach my $key (keys %sample_index) { 
			my @genotype = split ":", $arr[$key];
			if (/GT:AD:DP/) {
				next if (!$genotype[2] or $genotype[2] < 5);
			} else {
				next if (!$genotype[2] or $genotype[1] < 5);
			}
			$hash{$sample_index{$key}}{$genotype[0]}++;
		}
	}
	close VCF;
	foreach my $key (keys %sample_index) {
		my $het_ratio = sprintf("%.4f", $hash{$sample_index{$key}}{'0/1'}/($hash{$sample_index{$key}}{'0/0'} + $hash{$sample_index{$key}}{'1/1'}));
		my $sex;
		if ($het_ratio > 0.004) {
			$sex = 'Female';
		} else {
			$sex = 'Male';
		}
		print "$sample_index{$key}\t$het_ratio\t$sex\n";
	}
}
