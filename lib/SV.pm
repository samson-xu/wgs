package SV;

use File::Basename;
use WriteShell;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(manta);

sub manta {
	my $region = shift;
	my $manta = shift;
	my $convertInversion = shift;
	my $samtools = shift;
	my $AnnotSV = shift;
	my $svDB = shift;
	my $bam = shift;
	my $ref = shift;
	my $outDir = shift;
	my $thread = shift;
	my $bgzip = shift;
	my $tabix = shift;
	my $phenotypeDB = shift;
	my $iconv = shift;
	my $prefix = basename($bam);
	$prefix =~ s/.bam$//;
	$prefix =~ s/.final$//;
	my $manta_shell;
	my $target_flag = '';
	if ($region =~ m/.dict$/) {
		$manta_shell=<<RBED;

grep "^\@SQ" $region | cut -f 2,3 -d ":" | sed 's/LN://' | awk '{if(length(\$1) < 6) print \$1"\\t0\\t"\$2}' > $outDir/region.bed
RBED
	} else {
		$manta_shell="cp $region $outDir/region.bed\n";
		$target_flag = '--exome';

	}
	$region = "$outDir/region.bed";
	$manta_shell.=<<MANTA;
$bgzip -f $region
$tabix -f $region.gz
# manta shell
$manta \\
--bam $bam \\
--referenceFasta $ref \\
--runDir $outDir --callRegions $region.gz $target_flag 

$outDir/runWorkflow.py -j $thread

$convertInversion $samtools $ref $outDir/results/variants/diploidSV.vcf.gz > $outDir/$prefix.sv.vcf

$AnnotSV -annotationsDir $svDB -SVminSize 10 -SVinputFile $outDir/$prefix.sv.vcf -outputFile $outDir/$prefix.sv.annot

$phenotypeDB/sv_hpo.pl $outDir/$prefix.sv.annot.tsv $phenotypeDB/phenotype_hpo.txt $phenotypeDB/hpo_ch_info.txt > $outDir/$prefix.sv.annot.phenotype.tsv

$iconv -f utf-8 -t gb18030 $outDir/$prefix.sv.annot.phenotype.tsv > $outDir/$prefix.sv.annot.phenotype.ch.tsv

rm -rf $outDir/region* $outDir/workspace $outDir/run* 

MANTA

	write_shell($manta_shell, "$outDir/$prefix.manta.sh");
	
}
