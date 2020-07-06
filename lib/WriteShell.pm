package WriteShell;

use File::Basename;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(write_shell);

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

1;
