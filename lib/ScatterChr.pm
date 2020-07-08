package ScatterChr;

use POSIX;


require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(scatter_chr);

sub scatter_chr {
	my $file = shift;
	my $outDir = shift;
	my @intervals;
	if ($file =~ /dict$/) {
		my %chr_size;
		my @chrs;
		open D, $file or die $!;
		while (<D>) {
			next unless (/^\@SQ/);
			my @arr = split /\t/;
			$arr[1] =~ s/^SN://;
			$arr[2] =~ s/^LN://;
			$chr_size{$arr[1]} = $arr[2];
			push @chrs, $arr[1];
		}
		close D;
		my @sort_chr = sort {$chr_size{$b} <=> $chr_size{$a}} keys %chr_size;
		my $longest_sequence = $chr_size{$sort_chr[0]};
		my $tmp_size = 0;
		my $tmp_interval = '';
		for(my $i=0; $i<@chrs; $i++) {
			if ($chr_size{$chrs[$i]} + $tmp_size <= $longest_sequence) {
				$tmp_size += $chr_size{$chrs[$i]};
				$tmp_interval .= "-L $chrs[$i] ";
			} else {
				push @intervals, $tmp_interval;
				$tmp_interval = "-L $chrs[$i] ";
				$tmp_size = $chr_size{$chrs[$i]};
			}
			push @intervals, $tmp_interval if ($i == @chrs - 1);
		}
	} else {

	}
	return \@intervals;
}

1;
