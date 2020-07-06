package SampleStat;

use POSIX;


require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(stat_log);

sub stat_log {
	my $total = shift;
	my $dir = shift;
	my $time = strftime("%Y%m%d",localtime());
	system("mkdir -p $dir/stat_log") == 0 or die $!;
	open LOG, ">>$dir/stat_log/sample_stat.log" or die $!;
	print LOG "$time\t$total\n";
	close LOG;
}

1;
