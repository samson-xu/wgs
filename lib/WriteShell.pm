package WriteShell;

use File::Basename;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(write_shell parallel_shell);

sub write_shell {
	my $shell = shift;
	my $path = shift;
	my $dir = dirname($path);
	my $file = basename($path);
	system("mkdir -p $dir") == 0 || die $!;
	my $separator = "*" x 50;
	my $content=<<CONTENT;
#!/bin/bash
echo "${separator}$file Start at:`date '+%Y/%m/%d  %H:%M:%S'`$separator"
$shell	
echo "${separator}$file End at:`date '+%Y/%m/%d  %H:%M:%S'`$separator"
CONTENT
	open OUT, ">$path" or die $!;	
	print OUT $content;
	close OUT;
	system("chmod 755 $path") == 0 || die $!;
}

sub parallel_shell {
	my $shell = shift;
	my $path = shift;
	my $tt = shift;
	my $st = shift;
	my $max_parall_number = int($tt*1.2/$st);
	my $dir = dirname($path);
	my $file = basename($path);
	system("mkdir -p $dir") == 0 || die $!;
	my @shells = split /\n/, $shell;
	my ($new_shell, $i);
	my $j = 1;
	my $count = $max_parall_number;
	foreach my $s (@shells) {
		$i++;
		if ($i <= $count) {
			$new_shell .= "$s &\n";
		} else {
			$new_shell .= "wait\n";
			$new_shell .= "$s &\n";
			$j++;
			$count = $j * $max_parall_number;
		}
	}
	$new_shell .= "wait\n" unless ($new_shell =~ m/wait$/);
	my $separator = "*" x 50;
	my $content=<<CONTENT;
#!/bin/bash
echo "${separator}$file Start at:`date '+%Y/%m/%d  %H:%M:%S'`$separator"
$new_shell
echo "${separator}$file End at:`date '+%Y/%m/%d  %H:%M:%S'`$separator"
CONTENT
	open OUT, ">$path" or die $!;
	print OUT $content;
	close OUT;
	system("chmod 755 $path") == 0 || die $!;
}

1;
