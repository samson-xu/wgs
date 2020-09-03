package VariantCall;

use File::Basename;
use ScatterChr;
use WriteShell;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(call_variant);

sub call_variant {
	my $gatk = shift;
	my $bam = shift;
	my $ref = shift;
	my $dict = shift;
	my $dbsnp = shift;
	my $mills = shift;
	my $axiom = shift;
	my $hapmap = shift;
	my $omni = shift;
	my $thousand = shift;
	my $outDir = shift;
	my $target = shift;
	my $interval_padding = shift;
	if ($target) {
		$interval_padding = "-ip $interval_padding ";
	} else {
		$interval_padding = "";
	}
	my $prefix = basename($bam);
	$prefix =~ s/.final.bam$//;
	my ($gvcf_shell, $shell);	
	$gvcf_shell = "# HaplotypeCaller per-sample in GVCF mode\n";
	my @intervals;
	if ($target) {
		@intervals = @{scatter_normal_chr($target, "$outDir/split_bed")};
	} else {
		@intervals = @{scatter_normal_chr($dict)};
	}
	my ($groups_gvcf, $group_gvcf_file);
	foreach my $interval (@intervals) {
		my @arr = split /\s+/, $interval;
		my $group = basename($arr[1]);
		$groups_gvcf .= "--INPUT $outDir/$prefix.$group.g.vcf.gz \\\n";
		$group_gvcf_file .= "$outDir/$prefix.$group.g.vcf.gz* ";
		$gvcf_shell.=<<HaplotypeCaller;
$gatk HaplotypeCaller \\
-R $ref \\
-I $bam \\
$interval $interval_padding\\
-O $outDir/$prefix.$group.g.vcf.gz \\
-G AS_StandardAnnotation \\
-GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \\
-ERC GVCF &

HaplotypeCaller
	}

	$gvcf_shell.="wait\n\n";

	$gvcf_shell.=<<MergeGVCFs;
# Merge GVCFs generated per-interval for the same sample
$gatk MergeVcfs \\
$groups_gvcf--OUTPUT $outDir/$prefix.g.vcf.gz

rm $group_gvcf_file
MergeGVCFs
	write_shell($gvcf_shell, "$outDir/$prefix.gvcf.sh");	
	system("echo '$prefix\t$outDir/$prefix.g.vcf.gz' > $outDir/$prefix.gvcf.txt") == 0 || die $!;
	my $sample_gvcf_list_file = "$outDir/$prefix.gvcf.txt";
	my $output_vcf = genotype_gvcf($gatk, $sample_gvcf_list_file, \@intervals, $interval_padding, $outDir, $prefix, $ref, $dbsnp);
	vcf_filter($gatk, $output_vcf, $outDir, $prefix, '1', $dbsnp, $mills, $axiom, $hapmap, $omni, $thousand);
	$shell=<<SHELL;
# Pipeline for SNP/InDel detect
sh $outDir/$prefix.gvcf.sh >$outDir/$prefix.gvcf.sh.o 2>$outDir/$prefix.gvcf.sh.e
sh $outDir/$prefix.import.sh >$outDir/$prefix.import.sh.o 2>$outDir/$prefix.import.sh.e 
sh $outDir/$prefix.genotype.sh >$outDir/$prefix.genotype.sh.o 2>$outDir/$prefix.genotype.sh.e
sh $outDir/$prefix.vcf_filter.sh >$outDir/$prefix.vcf_filter.sh.o 2>$outDir/$prefix.vcf_filter.sh.e 
rm -rf $outDir/split_bed
rm $outDir/$prefix.vcf.gz
SHELL
	write_shell($shell, "$outDir/$prefix.variant.sh");

}

sub genotype_gvcf {
	my $gatk = shift;
	my $input = shift;
	my $intervals_ref = shift;
	my $interval_padding = shift;
	my $outDir = shift;
	my $prefix = shift;
	my $ref = shift;
	my $dbsnp = shift;
        my $import_shell = "# Import sample GVCFs\n";
        foreach my $interval (@{$intervals_ref}) {
                my @arr = split /\s+/, $interval;
                my $group = basename($arr[1]);
		$import_shell.=<<ImportGVCFs;		
$gatk --java-options -Xms8g \\
GenomicsDBImport \\
--genomicsdb-workspace-path $outDir/genomicsdb.$group \\
--batch-size 0 \\
$interval $interval_padding \\
--sample-name-map $input \\
--reader-threads 5 \\
--merge-input-intervals \\
--consolidate &

ImportGVCFs
	}	
	$import_shell .= "wait\n\n";
	write_shell($import_shell, "$outDir/$prefix.import.sh");
	my $genotype_shell = "# Genotype GVCFs\n";
        my ($groups_genotype, $group_genotype_file);
        foreach my $interval_genotype (@{$intervals_ref}) {
                my @arr = split /\s+/, $interval_genotype;
                my $group = basename($arr[1]);
                $groups_genotype .= "-I $outDir/$prefix.$group.vcf.gz \\\n";
                $group_genotype_file .= "$outDir/$prefix.$group.vcf.gz* ";
		$genotype_shell.=<<GenotypeGVCFs;
$gatk --java-options -Xms8g \\
GenotypeGVCFs \\
-R $ref \\
-O $outDir/$prefix.$group.vcf.gz \\
-D $dbsnp \\
-G AS_StandardAnnotation \\
--only-output-calls-starting-in-intervals \\
-V gendb://$outDir/genomicsdb.$group \\
$interval_genotype $interval_padding\\
--merge-input-intervals & 

GenotypeGVCFs
	}
	$genotype_shell .= "wait\n\n";
	$genotype_shell .= "rm -rf $outDir/genomicsdb*\n\n";
	$genotype_shell.=<<GatherVcfs;
# Gather split vcfs
$gatk --java-options -Xms6g \\
GatherVcfs \\
$groups_genotype-O $outDir/$prefix.vcf.gz

$gatk IndexFeatureFile \\
-I $outDir/$prefix.vcf.gz

rm $group_genotype_file

GatherVcfs
	write_shell($genotype_shell, "$outDir/$prefix.genotype.sh");
	return "$outDir/$prefix.vcf.gz" 
}

sub vcf_filter {
	my $gatk = shift;
	my $vcf = shift;
	my $outDir = shift;
	my $prefix = shift;
	my $flag = shift;
	my $dbsnp = shift;
	my $mills = shift;
	my $axiom = shift;
	my $hapmap = shift;
	my $omni = shift;
	my $thousand = shift;
	my $vcf_filter_shell = "# Filter VCF\n";
	if ($flag) {
		$vcf_filter_shell.=<<HardFilter;
# Hard Filter VCF
$gatk SelectVariants \\
-V $vcf \\
--select-type SNP \\
-O $outDir/$prefix.snp.vcf.gz & 

$gatk SelectVariants \\
-V $vcf \\
--select-type INDEL \\
-O $outDir/$prefix.indel.vcf.gz &

wait

$gatk VariantFiltration \\
-V $outDir/$prefix.snp.vcf.gz \\
--filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || QUAL < 30.0" \\
--filter-name "PASS" \\
--verbosity ERROR \\
-O $outDir/$prefix.snp.filter.vcf.gz &

$gatk VariantFiltration \\
-V $outDir/$prefix.indel.vcf.gz \\
--filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -20.0 || QUAL < 30.0" \\
--filter-name "PASS" \\
--verbosity ERROR \\
-O $outDir/$prefix.indel.filter.vcf.gz &

wait

$gatk MergeVcfs \\
-I $outDir/$prefix.snp.filter.vcf.gz \\
-I $outDir/$prefix.indel.filter.vcf.gz \\
-O $outDir/$prefix.filter.vcf.gz

$gatk IndexFeatureFile \\
-I $outDir/$prefix.filter.vcf.gz

rm $outDir/$prefix.snp.* $outDir/$prefix.indel.*

HardFilter
	} else {
		$vcf_filter_shell.=<<VQSR;
# VQSR Filter VCF
$gatk --java-options -Xms24g \\
VariantRecalibrator \\
-V $vcf \\
-O $outDir/$prefix.indel.recal.vcf.gz \\
--tranches-file $outDir/$prefix.indel.tranches \\
--trust-all-polymorphic \\
-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 \\
-tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 \\
-tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \\
-an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \\
--use-allele-specific-annotations \\
--mode INDEL \\
--max-gaussians 4 \\
-resource:mills,known=false,training=true,truth=true,prior=12 $mills \\
-resource:axiomPoly,known=false,training=true,truth=false,prior=10 $axiom \\
-resource:dbsnp,known=true,training=false,truth=false,prior=2 $dbsnp &

$gatk VariantRecalibrator \\
-V $vcf \\
-O $outDir/$prefix.snp.recal.vcf.gz \\
--tranches-file $outDir/$prefix.snp.tranches \\
--trust-all-polymorphic \\
-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 \\
-tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 \\
-tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \\
-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \\
--use-allele-specific-annotations \\
--mode SNP \\
--max-gaussians 6 \\
-resource:hapmap,known=false,training=true,truth=true,prior=15 $hapmap \\
-resource:omni,known=false,training=true,truth=true,prior=12 $omni \\
-resource:1000G,known=false,training=true,truth=false,prior=10 $thousand \\
-resource:dbsnp,known=true,training=false,truth=false,prior=7 $dbsnp &

wait

$gatk --java-options -Xms5g \\
ApplyVQSR \\
-O $outDir/$prefix.indel.filter.vcf.gz \\
-V $vcf \\
--recal-file $outDir/$prefix.indel.recal.vcf.gz \\
--use-allele-specific-annotations \\
--tranches-file $outDir/$prefix.indel.tranches \\
--truth-sensitivity-filter-level 99.0 \\
--create-output-variant-index true \\
--mode INDEL

$gatk --java-options -Xms5g \\
ApplyVQSR \\
-O $outDir/$prefix.filter.vcf.gz \\
-V $outDir/$prefix.indel.filter.vcf.gz \\
--recal-file $outDir/$prefix.snp.recal.vcf.gz \\
--use-allele-specific-annotations \\
--tranches-file $outDir/$prefix.snp.tranches \\
--truth-sensitivity-filter-level 99.7 \\
--create-output-variant-index true \\
--mode SNP

rm $outDir/$prefix.snp.* $outDir/$prefix.indel.*

VQSR
	}
	write_shell($vcf_filter_shell, "$outDir/$prefix.vcf_filter.sh");
}

1;
