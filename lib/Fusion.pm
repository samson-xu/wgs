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
	my $region = shift;
	my $arg = shift;
	my $sampleId = shift;
	if ($region =~ m/.dict$/) {
	$factera_shell=<<RBED;

grep "^\@SQ" $region | cut -f 2,3 -d ":" | sed 's/LN://' | awk '{if(length(\$1) < 6) print \$1"\\t0\\t"\$2}' > $outDir/region.bed
RBED
        } else {
                $factera_shell="cp $region $outDir/region.bed\n";
        }
        $region = "$outDir/region.bed";

	my $factera_shell.=<<FUSION;
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
       
mv $outDir/$sampleId.final.factera.fusions.txt $outDir/$sampleId.final.factera.fusions.xls
       
rm $outDir/*blast*
FUSION

	write_shell($factera_shell, "$outDir/$sampleId.fusion.sh");

}