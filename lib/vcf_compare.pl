#!/usr/bin/env perl
#use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;
use Getopt::Long;
use Data::Dumper;

# Global variable
my ($help, $outDir, %site, $file_number, $var_same, $snp_same, $indel_same, %store);

($var_same, $snp_same, $indel_same) = (0, 0, 0);

# Get Parameter
GetOptions(
	"h|help" => \$help,	
	"outDir=s" => \$outDir,
);

my $tmpDir = `pwd`;
chomp $tmpDir;
$outDir ||= "$tmpDir";

# Guide for program
my $guide_separator = "#" x 80;
my $program = basename(abs_path($0));
$program =~ s/\.pl//;
my $guide=<<INFO;
VER
	AUTHOR: xuxiangyang(xy_xu\@foxmail.com)
	NAME: $program
	PATH: $0
        VERSION: v0.1	2016-10-26
NOTE
	1. Input file should be more than 2

USAGE
	$program <options> var1.file,var2.file... [varfile.lst]
	$guide_separator Basic $guide_separator
	--help			print help information
	--outDir <s>		script out Dir, default "$outDir"

INFO
# Call guide
die $guide if (@ARGV == 0 || defined $help);

# Program check

# Main
my $input = shift;
if ($input =~ /,/i) {
	my @files = split ",", $input;
	$file_number = @files;
	foreach my $file (@files) {
		my ($var_number, $snp_number, $indel_number);
		($var_number, $snp_number, $indel_number) = (0, 0, 0);
		my $dir = `dirname $file`;
		chomp $dir;
		my $file = basename($file);
		-s "$dir/$file" or die $!;
		if ($file =~ /gz$/) {
			open FL, "gzip -dc $dir/$file |" or die $!;
		} else {
			open FL, "$dir/$file" or die $!;
		}
		while (<FL>) {
			next if (/^#/);
			$var_number++;
			chomp;
			my @row = split /\t/;
			if (length($row[3]) ne length($row[4])) {
				$indel_number++;	
			} else {
				$snp_number++;
			}
			if ($site{$row[0]}{$row[1]}{$row[3]}{$row[4]}) {
				$site{$row[0]}{$row[1]}{$row[3]}{$row[4]}++;
			} else {
				$site{$row[0]}{$row[1]}{$row[3]}{$row[4]} = 1;
			}	
		}
		close FL;
		push @{$store{"$file"}}, ($var_number, $snp_number, $indel_number);
	}
} else {
	unless ($input =~ m/\.lst/i) {
		print "$input may be not a list!\n";
		exit;
	}	
	open INPUT, $input or die $!;
	while (<INPUT>)	 {
		next if (/^#/);
		$file_number++;
		chomp;
		$file = $_;
		my $dir = `dirname $file`;
		chomp $dir;
                my $file = basename($file);
		my ($var_number, $snp_number, $indel_number);
		($var_number, $snp_number, $indel_number) = (0, 0, 0);
		-s "$dir/$file" or die $!;
		if ($file =~ /gz$/) {
            open FL, "gzip -dc $dir/$file |" or die $!;
		} else {
			open FL, "$dir/$file" or die $!;
		}
        while (<FL>) {
				next if (/^#/);
                $var_number++;
                chomp;
                my @row = split /\t/;
                if (length($row[3]) ne length($row[4])) {
                        $indel_number++;
                } else {
                        $snp_number++;
                }
                if ($site{$row[0]}{$row[1]}{$row[3]}{$row[4]}) {
                        $site{$row[0]}{$row[1]}{$row[3]}{$row[4]}++;
                } else {
                        $site{$row[0]}{$row[1]}{$row[3]}{$row[4]} = 1;
                }
		}
		close FL;
		push @{$store{"$file"}}, ($var_number, $snp_number, $indel_number);
	}
	close INPUT;
}

foreach my $chr (keys %site) {
	foreach my $start (keys %{$site{$chr}}) {
		foreach my $ref (keys %{$site{$chr}{$start}}) {
			foreach my $alt (keys %{$site{$chr}{$start}{$ref}}) {	
				if ($site{$chr}{$start}{$ref}{$alt} == $file_number) {
					$var_same++;
					if (length($ref) ne length($alt)) {
						$indel_same++;
					} else {
						$snp_same++;
					}
				}					
			}
		}
	}
}

printf "%-20s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n", "File_name", "Var_number", "Var_same", "Var_rate(%)", "SNP_number", "SNP_same", "SNP_rate(%)", "InDel_number", "InDel_same", "InDel_rate(%)";

foreach my $var_file (sort {$a cmp $b} keys %store) {
	if ($store{$var_file}[0] == 0) {
		print STDERR "There is no variation in the $var_file!\n";
		next;
	}
	if ($store{$var_file}[1] == 0) {
		print STDERR "There is no SNP in the $var_file!\n";
		printf "%-20s\t%10d\t%10d\t%10.2f\t%10d\t%10d\t%10.2f\t%10d\t%10d\t%10.2f\n", $var_file, $store{$var_file}[0], $var_same, $var_same/$store{$var_file}[0]*100, $store{$var_file}[1], $snp_same, 0, $store{$var_file}[2], $indel_same, $indel_same/$store{$var_file}[2]*100;
	}
	if ($store{$var_file}[2] == 0) {
		print STDERR "There is no InDel in the $var_file!\n";
		printf "%-20s\t%10d\t%10d\t%10.2f\t%10d\t%10d\t%10.2f\t%10d\t%10d\t%10.2f\n", $var_file, $store{$var_file}[0], $var_same, $var_same/$store{$var_file}[0]*100, $store{$var_file}[1], $snp_same, $snp_same/$store{$var_file}[1]*100, $store{$var_file}[2], $indel_same, 0;
	}
	if ($store{$var_file}[0] != 0 && $store{$var_file}[1] != 0 && $store{$var_file}[2] != 0) {
		printf "%-20s\t%10d\t%10d\t%10.2f\t%10d\t%10d\t%10.2f\t%10d\t%10d\t%10.2f\n", $var_file, $store{$var_file}[0], $var_same, $var_same/$store{$var_file}[0]*100, $store{$var_file}[1], $snp_same, $snp_same/$store{$var_file}[1]*100, $store{$var_file}[2], $indel_same, $indel_same/$store{$var_file}[2]*100;
	}
}
