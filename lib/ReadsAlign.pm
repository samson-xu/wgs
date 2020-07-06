package ReadsAlign;

use File::Basename;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(reads_align);

sub reads_align {
	my $program = shift;
	my $samtools = shift;
	my $gatk = shift; 
	my $thread = shift;
	my $fastq = shift;
	my $ref = shift;
	my $way = shift;
	my $arg = shift;
	my $outDir = shift;
	my $prefix = basename($fastq->[0]);
	$prefix =~ s/.clean.1.fq.gz//;
	my $shell=<<FastqToSam;
# Convert FASTQ to read unmapped BAM
$gatk FastqToSam \\
-F1 $fastq->[0] \\
-F2 $fastq->[1] \\
-O $outDir/$prefix.unmapped.bam \\
-SM $prefix \\
-PL ILLUMINA \\
-RG $prefix &

# BWA for alignment
FastqToSam
	if ($way eq 'backtrack') {
		$arg = "-t $thread" if ($arg eq '');
		$shell .= "$program aln $arg $ref $fastq->[0] > $outDir/$prefix.aln_sa1.sai &\n";
		$shell .= "$program aln $arg $ref $fastq->[1] > $outDir/$prefix.aln_sa2.sai &\n";
		$shell .= "wait\n";
		$shell .= "$program sampe $ref $outDir/$prefix.aln_sa1.sai $outDir/$prefix.aln_sa2.sai $fastq->[0] $fastq->[1] | $samtools view -1 --threads $thread - > $outDir/$prefix.bam\n\n";
	} else {
		$arg = "-K 100000000 -t $thread -Y" if ($arg eq '');
		$shell .= "$program mem $arg $ref $fastq->[0] $fastq->[1] | $samtools view -1 --threads $thread - > $outDir/$prefix.bam &\n";
		$shell .= "wait\n\n";
	}
	$shell.=<<MergeBam;
# Merge original input uBAM file with BWA-aligned BAM file
$gatk MergeBamAlignment \\
--VALIDATION_STRINGENCY SILENT \\
--EXPECTED_ORIENTATIONS FR \\
--ATTRIBUTES_TO_RETAIN X0 \\
--ALIGNED_BAM $outDir/$prefix.bam \\
--UNMAPPED_BAM $outDir/$prefix.unmapped.bam \\
--OUTPUT $outDir/$prefix.merge.bam \\
--REFERENCE_SEQUENCE $ref \\
--SORT_ORDER "unsorted" \\
--IS_BISULFITE_SEQUENCE false \\
--ALIGNED_READS_ONLY false \\
--CLIP_ADAPTERS false \\
--MAX_RECORDS_IN_RAM 2000000 \\
--ADD_MATE_CIGAR true \\
--MAX_INSERTIONS_OR_DELETIONS -1 \\
--PRIMARY_ALIGNMENT_STRATEGY MostDistant \\
--UNMAPPED_READ_STRATEGY COPY_TO_TAG \\
--ALIGNER_PROPER_PAIR_FLAGS true \\
--UNMAP_CONTAMINANT_READS true

MergeBam

	$shell.=<<MarkDuplicates;
# Mark duplicate reads to avoid counting non-independent observations
$gatk MarkDuplicates \\
--INPUT $outDir/$prefix.merge.bam \\
--OUTPUT $outDir/$prefix.markdup.bam \\
--METRICS_FILE $outDir/$prefix.metrics \\
--VALIDATION_STRINGENCY SILENT \\
--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \\
--ASSUME_SORT_ORDER "queryname" \\
--CREATE_MD5_FILE false

MarkDuplicates

	$shell.=<<SortSam;
# Sort BAM file by coordinate order and fix tag values for NM and UQ
$gatk SortSam \\
--INPUT $outDir/$prefix.markdup.bam \\
--OUTPUT /dev/stdout \\
--SORT_ORDER "coordinate" \\
--CREATE_INDEX false \\
--CREATE_MD5_FILE false \\
| \\
$gatk SetNmMdAndUqTags \\
--INPUT /dev/stdin \\
--OUTPUT $outDir/$prefix.markdup.sort.bam \\
--CREATE_INDEX true \\
--CREATE_MD5_FILE true \\
--REFERENCE_SEQUENCE $ref

SortSam
	

	return $shell;
}

1;
