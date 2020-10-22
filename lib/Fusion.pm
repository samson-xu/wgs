package Fusion;

use File::Basename;
use WriteShell;

#write_shell($manta_shell, "$outDir/$prefix.manta.sh");

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(factera);

sub factera {
	my $fusion = shift;
	my $outDir = shift;
	my $samtools = shift;	
	my $twoBitToFa = shift;
	my $blastn = shift;
	my $makeblastdb = shift;
	my $thread = shift;
	my $bam = shift;
	my $interGene = shift;
	my $bit = shift;
	my $phenotype = shift;
	my $iconv = shift;
	my $region = shift;
	my $arg = shift;
	my $sampleId = shift;
	my $factera_shell;
	if ($region =~ m/.dict$/) {
		$factera_shell=<<RBED;
grep "^\@SQ" $region | cut -f 2,3 -d ":" | sed 's/LN://' | awk '{if(length(\$1) < 6) print \$1"\\t0\\t"\$2}' > $outDir/region.bed
RBED
        } else {
                $factera_shell = "cp $region $outDir/region.bed\n";
        }
        $region = "$outDir/region.bed";
	$factera_shell.=<<FUSION;
# Fusion detect
$fusion \\
-o $outDir \\
-A $samtools \\
-T $twoBitToFa \\
-B $blastn \\
-M $makeblastdb \\
-p $thread \\
-F $arg \\
$bam \\
$interGene \\
$bit $region
       
$phenotype/fusion_hpo.pl $outDir/$sampleId.final.factera.fusions.txt $phenotype/omim_gene_phenotypes.txt $phenotype/phenotype_hpo.txt $phenotype/hpo_ch_info.txt > $outDir/$sampleId.final.factera.fusions.xls

$iconv -f utf-8 -t gb18030 $outDir/$sampleId.final.factera.fusions.xls > $outDir/$sampleId.final.factera.fusions.ch.xls
       
rm $outDir/*blast* $outDir/region.bed $outDir/*discordantpair*
FUSION

	write_shell($factera_shell, "$outDir/$sampleId.fusion.sh");

}
