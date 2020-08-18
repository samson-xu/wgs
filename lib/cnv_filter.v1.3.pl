#!/usr/bin/env perl
use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

# Global variable
my ($help, $outDir, $output_file);

# Get Parameter
GetOptions(
    "h|help" => \$help, 
    "outDir=s" => \$outDir,
    "output_file=s" => \$output_file
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
    AUTHOR: samson-xu(xy_xu\@foxmail.com)
    NAME: $program
    PATH: $0
    VERSION: v1.3   2020-08-10
NOTE
    Filter CNV according to sex information

USAGE
    $program <options> sex.txt *.cnv 
    $guide_separator Basic $guide_separator
    --help              print help information
    --outDir <s>        script out Dir, default "$outDir"
    --output_file <s>	output file name	

INFO

die $guide if (@ARGV == 0 || defined $help);

# Main
my $prefix = basename($ARGV[1]);
$prefix =~ s/.cns.call//;
#print "$prefix\n";
my $sex;

open SEX, $ARGV[0] or die $!;
open CNV, $ARGV[1] or die $!;

while (<SEX>) {
    chomp;
    next if (/^\s*$/);
    next if (/^#/);
    if (/$prefix.cnr/) {
        my @arr = split;
        $sex = $arr[1];
        last;
    }	
}
close SEX;

while (<CNV>) {
    chomp;
    next if (/^\s*$/);
    next if (/^#/);
    next if (/^chromosome/);
    my $cnv_status;
    my @arr = split;
    my $cn = $arr[7];
    if ($arr[0] =~ m/chrx/i) {
        if ($sex eq "Male") {
            next if ($cn == 1);
            if ($cn > 1) {
                $cnv_status = 'DUP';
            } else {
                $cnv_status = 'DEL';
            }
        } else {
            next if ($cn == 2);
            if ($cn > 2) {
                $cnv_status = 'DUP';
            } else {
                $cnv_status = 'DEL';
            }
        }
    } elsif ($arr[0] =~ m/chry/i) {
        if ($sex eq "Male") {
            next if ($cn == 1);
            if ($cn > 1) {
                $cnv_status = 'DUP';
            } else {
                $cnv_status = 'DEL';
            }
        } else {
            next;
        }
    } else {
        next if ($cn == 2);
        if ($cn > 2) {
            $cnv_status = 'DUP';
        } else {
            $cnv_status = 'DEL';
        }
    }
    if ($output_file) {
        open OUT, ">$output_file" or die $!;
        print OUT "$arr[0]\t$arr[1]\t$arr[2]\t$cn\t$cnv_status\n";
        close OUT;
    } else {
        print "$arr[0]\t$arr[1]\t$arr[2]\t$cn\t$cnv_status\n";
    }
}
close CNV;

