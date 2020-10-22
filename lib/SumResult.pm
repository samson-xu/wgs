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
	$copy_shell .= "cp $inDir/03.snp_indel/$sampleId.*vcf.gz $outDir\n" if ($step =~ m/v/);
	$copy_shell .= "cp $inDir/04.sv/manta/*.phenotype.ch.tsv $outDir\n" if ($step =~ m/s/);
	$copy_shell .= "cp $inDir/04.sv/fusion/*.ch.xls $outDir\n" if ($step =~ m/u/);
	$copy_shell .= "cp $inDir/05.mtdna/*.xls $outDir\n" if ($step =~ m/m/);
	#$copy_shell .= "cp $inDir/../*cnv/cnv/$sampleId*tsv $outDir\n" if ($step =~ m/n/);
	write_shell($copy_shell, "$inDir/$sampleId.summary.sh");
}
