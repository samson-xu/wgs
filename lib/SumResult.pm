package SumResult;

use File::Basename;
use WriteShell;

#write_shell($manta_shell, "$outDir/$prefix.manta.sh");

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(copy);

sub copy {
	my $inDir = shift;
	my $outDir = shift;
	my $step = shift;
	my $sampleId = shift;
	system("mkdir -p $outDir") == 0 || die $!;
	my $copy_shell = "#script for copy result file\n";
	my $copy_shell .= "cp $inDir/03.snp_indel/$sampleId.*.vcf.gz $outDir" if ($step =~ m/v/);
	my $copy_shell .= "cp $inDir/04.sv/manta/*.tsv $outDir" if ($step =~ m/s/);
	my $copy_shell .= "cp $inDir/04.sv/fusion/*.xls $outDir" if ($step =~ m/u/);
	my $copy_shell .= "cp $inDir/05.mtdna/*.xls $outDir" if ($step =~ m/m/);
	my $copy_shell .= "cp $inDir/../*cnv/cnv/$sampleId*tsv $outDir" if ($step =~ m/n/);
	write_shell($copy_shell, "$inDir/$sampleId.summary.sh");
}
