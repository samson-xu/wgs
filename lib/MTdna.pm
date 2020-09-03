package MTdna;

use File::Basename;
use WriteShell;

#write_shell($manta_shell, "$outDir/$prefix.manta.sh");

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(mtdna);

sub mtdna {
	my $bam = shift;	
	my $gatk = shift;
	my $bwa = shift;
	my $thread = shift;
	my $annotDir = shift;
	my $binDir = shift;
	my $dbDir = shift;
	my $outDir = shift;
	my $convert2annovar = shift;
	my $table_annovar = shift;
	my $prefix = basename($bam);
	$prefix =~ s/.bam$//;
	$prefix =~ s/.final$//;
	my $mt_align = align("$outDir/$prefix.chrM.unmapped.bam", $gatk, $bwa, $outDir, $thread, "$dbDir/Homo_sapiens_assembly38.chrM.fasta", "normal");
	my $shift_align = align("$outDir/$prefix.chrM.unmapped.bam", $gatk, $bwa, $outDir, $thread, "$dbDir/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta", "shift");
	my $mt_call = mutect2("$outDir/$prefix.chrM.normal.bam", $gatk, "$dbDir/Homo_sapiens_assembly38.chrM.fasta", $outDir, "-L chrM:576-16024");
	my $shift_call = mutect2("$outDir/$prefix.chrM.shift.bam", $gatk, "$dbDir/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta", $outDir, "-L chrM:8025-9144");
	my $mtdna_shell=<<MT;
# Workflow for SNP and INDEL variant calling on mitochondria
# Subsets a whole genome bam to just Mitochondria reads
$gatk PrintReads \\
-L chrM \\
--read-filter MateOnSameContigOrNoMappedMateReadFilter \\
--read-filter MateUnmappedAndUnmappedReadFilter \\
-I $bam \\
-O $outDir/$prefix.subset.bam

# Removes alignment information while retaining recalibrated base qualities and original alignment tags
$gatk --java-options -Xmx1000m \\
RevertSam \\
--INPUT $outDir/$prefix.subset.bam \\
--OUTPUT_BY_READGROUP false \\
--OUTPUT $outDir/$prefix.chrM.unmapped.bam \\
--VALIDATION_STRINGENCY LENIENT \\
--ATTRIBUTE_TO_CLEAR FT \\
--ATTRIBUTE_TO_CLEAR CO \\
--SORT_ORDER queryname \\
--RESTORE_ORIGINAL_QUALITIES false

rm $outDir/$prefix.subset.*

$mt_align

$shift_align

wait

rm $outDir/$prefix.chrM.unmapped.bam  $outDir/$prefix.chrM.*.merge.bam $outDir/$prefix.chrM.*.mark.bam *.metrics 

# Collect coverage metrics
#$gatk --java-options -Xms2000m \\
#CollectWgsMetrics \\
#--INPUT $outDir/$prefix.chrM.normal.bam \\
#--VALIDATION_STRINGENCY SILENT \\
#--REFERENCE_SEQUENCE $dbDir/Homo_sapiens_assembly38.chrM.fasta \\
#--OUTPUT $outDir/metrics.txt \\
#--USE_FAST_ALGORITHM true \\
#--READ_LENGTH 151 \\
#--COVERAGE_CAP 100000 \\
#--INCLUDE_BQ_HISTOGRAM true \\
#--THEORETICAL_SENSITIVITY_OUTPUT $outDir/theoretical_sensitivity.txt
#
#R --vanilla <<CODE
#df = read.table("$outDir/metrics.txt",skip=6,header=TRUE,stringsAsFactors=FALSE,sep='\t',nrows=1)
#write.table(floor(df[,"MEAN_COVERAGE"]), "$outDir/mean_coverage.txt", quote=F, col.names=F, row.names=F)
#CODE

#Uses Haplochecker to estimate levels of contamination in mitochondria
java -jar $binDir/mitolib-0.1.2.jar haplochecker \\
--in $outDir/$prefix.chrM.normal.bam \\
--ref $dbDir/Homo_sapiens_assembly38.chrM.fasta \\
--out $outDir/haplochecker_out \\
--QUAL 20 \\
--MAPQ 30 \\
--VAF 0.01

#python3 <<CODE
#
#import csv
#
#with open("$outDir/haplochecker_out/$prefix.chrM.normal.contamination.txt") as output:
#    reader = csv.DictReader(output, delimiter='\t')
#    for row in reader:
#        print(row["MajorHG"], file=open("$outDir/major_hg.txt", 'w'))
#        print(row["MajorLevel"], file=open("$outDir/major_level.txt", 'w'))
#        print(row["MinorHG"], file=open("$outDir/minor_hg.txt", 'w'))
#        print(row["MinorLevel"], file=open("$outDir/minor_level.txt", 'w'))
#CODE

$mt_call

$shift_call 

wait

rm $outDir/*.log

# Lifts over shifted vcf of control region and combines it with the rest of the chrM calls
$gatk LiftoverVcf \\
-I $outDir/$prefix.chrM.shift.vcf \\
-O $outDir/$prefix.chrM.shifted_back.vcf \\
-R $dbDir/Homo_sapiens_assembly38.chrM.fasta \\
--CHAIN $dbDir/ShiftBack.chain \\
--REJECT $outDir/$prefix.chrM.rejected.vcf

$gatk MergeVcfs \\
-I $outDir/$prefix.chrM.shifted_back.vcf \\
-I $outDir/$prefix.chrM.normal.vcf \\
-O $outDir/$prefix.chrM.merge.vcf

# MergeStats
$gatk MergeMutectStats \\
--stats $outDir/$prefix.chrM.shift.vcf.stats \\
--stats $outDir/$prefix.chrM.normal.vcf.stats \\
-O $outDir/$prefix.chrM.raw.combined.stats

# Mutect2 Filtering for calling Snps and Indels
contamination=`awk '{print \$8}' $outDir/haplochecker_out/$prefix.chrM.normal.contamination.txt | grep -v MinorLevel`
$gatk --java-options -Xmx2500m \\
FilterMutectCalls \\
-V $outDir/$prefix.chrM.merge.vcf \\
-R $dbDir/Homo_sapiens_assembly38.chrM.fasta \\
-O $outDir/$prefix.chrM.merge.filtered.vcf \\
--stats $outDir/$prefix.chrM.raw.combined.stats \\
--max-alt-allele-count 4 \\
--mitochondria-mode true \\
--contamination-estimate \$contamination

$gatk VariantFiltration -V $outDir/$prefix.chrM.merge.filtered.vcf \\
-O $outDir/$prefix.chrM.final.vcf \\
--mask $dbDir/blacklist_sites.hg38.chrM.bed \\
--mask-name "blacklisted_site"

rm $outDir/$prefix.chrM.merge.* $outDir/$prefix.chrM.shift* $outDir/$prefix.chrM.normal*vcf* $outDir/*.stats
rm -rf haplochecker_out/

# Annotation for SNP and INDEL variant on mitochondria
grep -v "#" $outDir/$prefix.chrM.final.vcf | grep -v -E 'blacklisted_site|weak_evidence' > $outDir/$prefix.chrM.final.filter.vcf
$convert2annovar --format vcf4 $outDir/$prefix.chrM.final.filter.vcf > $outDir/$prefix.chrM.final.filter.av
$table_annovar --buildver hg19 --remove --protocol ensGene,mitomap20190903,mitimpact3 --operation g,f,f --nastring . $outDir/$prefix.chrM.final.filter.av $annotDir --outfile $outDir/$prefix.chrM.final
$binDir/add_AF.pl $outDir/$prefix.chrM.final.filter.vcf $outDir/$prefix.chrM.final.hg19_multianno.txt > $outDir/$prefix.chrM.final.result.xls
rm $outDir/$prefix.chrM.final.filter.vcf $outDir/$prefix.chrM.final.filter.av $outDir/$prefix.chrM.final.hg19_multianno.txt
rm -rf $outDir/haplochecker_out
MT

	write_shell($mtdna_shell, "$outDir/$prefix.mtdna.sh");

}

sub align {
	my $bam = shift;
	my $gatk = shift;
	my $bwa = shift;
	my $outDir = shift;
	my $thread = shift;
	my $ref = shift;
	my $flag = shift;
	my $prefix = basename($bam);
	$prefix =~ s/.unmapped.bam$//;
	$prefix .= ".$flag";
	my $align_shell=<<ALIGN;
# Aligns with BWA and MergeBamAlignment, then Marks Duplicates. Outputs a coordinate sorted bam.
$gatk --java-options -Xms5000m \\
SamToFastq \\
--INPUT $bam \\
--FASTQ /dev/stdout \\
--INTERLEAVE true \\
--NON_PF true | \\
$bwa mem -K 100000000 -p -v 3 -t $thread -Y $ref /dev/stdin - 2> >(tee $outDir/$prefix.bwa.stderr.log >&2) | \\
$gatk --java-options -Xms3000m \\
MergeBamAlignment \\
--VALIDATION_STRINGENCY SILENT \\
--EXPECTED_ORIENTATIONS FR \\
--ATTRIBUTES_TO_RETAIN X0 \\
--ATTRIBUTES_TO_REMOVE NM \\
--ATTRIBUTES_TO_REMOVE MD \\
--ALIGNED_BAM /dev/stdin \\
--UNMAPPED_BAM $bam \\
--OUTPUT $outDir/$prefix.merge.bam \\
--REFERENCE_SEQUENCE $ref \\
--PAIRED_RUN true \\
--SORT_ORDER unsorted \\
--IS_BISULFITE_SEQUENCE false \\
--ALIGNED_READS_ONLY false \\
--CLIP_ADAPTERS false \\
--MAX_RECORDS_IN_RAM 2000000 \\
--ADD_MATE_CIGAR true \\
--MAX_INSERTIONS_OR_DELETIONS -1 \\
--PRIMARY_ALIGNMENT_STRATEGY MostDistant \\
--UNMAPPED_READ_STRATEGY COPY_TO_TAG \\
--ALIGNER_PROPER_PAIR_FLAGS true \\
--UNMAP_CONTAMINANT_READS true \\
--ADD_PG_TAG_TO_READS false \\
&& \\
$gatk --java-options -Xms4000m \\
MarkDuplicates \\
--INPUT $outDir/$prefix.merge.bam \\
--OUTPUT $outDir/$prefix.mark.bam \\
--METRICS_FILE $outDir/$prefix.metrics \\
--VALIDATION_STRINGENCY SILENT \\
--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \\
--ASSUME_SORT_ORDER queryname \\
--CLEAR_DT false \\
--ADD_PG_TAG_TO_READS false \\
&& \\
$gatk --java-options -Xms4000m \\
SortSam \\
--INPUT $outDir/$prefix.mark.bam \\
--OUTPUT $outDir/$prefix.bam \\
--SORT_ORDER coordinate \\
--CREATE_INDEX true \\
--MAX_RECORDS_IN_RAM 300000 &
ALIGN
	return $align_shell;
}

sub mutect2 {
	my $bam = shift;
	my $gatk = shift;
	my $ref = shift;
	my $outDir = shift;
	my $region = shift;
	my $prefix = basename($bam);
	$prefix =~ s/.bam$//;
	my $mutect2_shell=<<MUTECT;
# Mutect2 for calling Snps and Indels
$gatk Mutect2 \\
-R $ref \\
-I $bam \\
-O $outDir/$prefix.vcf \\
$region \\
--annotation StrandBiasBySample \\
--mitochondria-mode true \\
--max-reads-per-alignment-start 75 \\
--max-mnp-distance 0 &
MUTECT
	return $mutect2_shell;

}
