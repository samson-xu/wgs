package ScatterChr;

use POSIX;


require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(scatter_chr scatter_normal_chr);

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
		my %chr_size;
		my @chrs = map {"chr$_"} (1..22, 'X', 'Y');
		my %chr_region;
		open F, $file or die $!;
		while (<F>) {
			next if (/#/);	
			next if (/^\s*$/);
			chomp;
			my @arr = split /\s+/;
			$chr_size{$arr[0]} += $arr[2] - $arr[1];
			$chr_region{$arr[0]} .= "$_\n";
		}
		close F;
		system("mkdir -p $outDir") == 0 || die $!;
		foreach my $chr (keys %chr_region) {
			open OUT, ">$outDir/$chr.bed" or die $!; 
			print OUT $chr_region{$chr};
			close OUT;
		}
		my @sort_chr = sort {$chr_size{$b} <=> $chr_size{$a}} keys %chr_size;
		my $longest_sequence = $chr_size{$sort_chr[0]};
		my $tmp_size = 0;
		my $tmp_interval = '';
		for(my $i=0; $i<@chrs; $i++) {
			if ($chr_size{$chrs[$i]} + $tmp_size <= $longest_sequence) {
				$tmp_size += $chr_size{$chrs[$i]};
				$tmp_interval .= "-L $outDir/$chrs[$i].bed ";
			} else {
				push @intervals, $tmp_interval;
				$tmp_interval = "-L $outDir/$chrs[$i].bed ";
				$tmp_size = $chr_size{$chrs[$i]};
			}
			push @intervals, $tmp_interval if ($i == @chrs - 1);
		}
	}
	return \@intervals;
}

sub scatter_normal_chr {
	my $file = shift;
	my $outDir = shift;
	my @intervals;
	if ($file =~ /dict$/) {
		my @chrs;
		open D, $file or die $!;
		while (<D>) {
			next unless (/^\@SQ/);
			my @arr = split /\t/;
			$arr[1] =~ s/^SN://;
			$arr[2] =~ s/^LN://;
			next if ($arr[2] < 20000000);
			push @chrs, $arr[1];
		}
		close D;
		foreach my $chr (@chrs) {
			push @intervals, "-L $chr";
		}
	} else {
		my @chrs = map {"chr$_"} (1..22, 'X', 'Y');
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
			push @intervals, "-L $outDir/$chr.bed" if (-e "$outDir/$chr.bed");
		}
	}
	return \@intervals;
}

1;
