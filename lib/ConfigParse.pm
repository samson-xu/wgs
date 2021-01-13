package ConfigParse;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(path_check);

sub path_check {
	my $config_file = shift;
	my %config = ();
	open CON, $config_file or die $!;
	while (<CON>) {
		next if (/#/);		
		next if (/^\s*$/);
		chomp;
		my @arr = split /\s+/;
		if ($arr[1] and -e $arr[1]) {
			$config{$arr[0]} = $arr[1];
		} else {
			print STDERR "ERROR: $arr[0] does not exist or path error!\n";
		}
	}
	close CON;
	return \%config;
}

1;
