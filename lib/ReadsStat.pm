package ReadsStat;

use File::Basename;
use JSON;
use Data::Dumper;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(reads_stat);

# Command line function call example
# perl -I lib -MReadsStat -e "reads_stat('*.fastq.json')"

sub reads_stat {
	my $file = shift;
	my $outDir = dirname($file);
	my $prefix = basename($file);
	$prefix =~ s/.fastq.json//;
	my $context;
	open JSON, $file or die $!;
	while (<JSON>) {
	    $context .= $_;
	}
	close JSON;
	my $input_json = decode_json($context);
	my $output=<<OUTPUT; 
Sample\t$prefix
Clean_reads\t$input_json->{'summary'}->{'after_filtering'}->{'total_reads'}
Clean_bases\t$input_json->{'summary'}->{'after_filtering'}->{'total_bases'}
Q20_rate\t$input_json->{'summary'}->{'after_filtering'}->{'q20_rate'}
Q30_rate\t$input_json->{'summary'}->{'after_filtering'}->{'q30_rate'}
GC_content\t$input_json->{'summary'}->{'after_filtering'}->{'gc_content'}
Duplication_rate\t$input_json->{'duplication'}->{'rate'}
Insert_size\t$input_json->{'insert_size'}->{'peak'}
OUTPUT
	open OUT, ">$outDir/$prefix.fq.stat.txt" or die $!;
	print OUT $output;
	close OUT;
}

1;
