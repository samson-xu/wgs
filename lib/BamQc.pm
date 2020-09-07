package BamQc;

use POSIX;
use Parallel::ForkManager;
use File::Basename;
use Data::Dumper;


require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(bam_qc);

# Command line function call example
# perl -I lib -MBamQc -e "bam_qc('bam', 'samtools')"

sub bam_qc {
	my $bam = shift;
	my $samtools = shift;
	my $target = shift;
	my $bedtools = shift;
	my $thread = shift;
	my $interval_padding;
	$thread ||= 35;	
	$interval_padding ||= 200;
	my $outDir = dirname($bam);
	my $sampleId = basename($bam);
	$sampleId =~ s/.bam//;
	$sampleId =~ s/.final//;
	my (%depth_info, %bam_info, @all_chr, $output);
	my $bam_head = `$samtools view -H $bam | grep '^\@SQ' | awk '{print \$2}' | cut -f 2 -d ":"`;
	chomp $bam_head;
	if ($bam_head) {
		@all_chr = split(/\n/, $bam_head);
	} else {
		@all_chr = map {"chr$_"} (1..22, 'X', 'Y');		
	}
	my @all_chr_with_unmapped = ('unmapped', @all_chr);
	my $pm = new Parallel::ForkManager($thread);
	#%bam_info = %{bam_stat($bam, $samtools)};
	print "Start mapping stat by chr!\n";
	system(">$outDir/$sampleId.chr.stat.txt") == 0 || die $!;
	foreach my $part (@all_chr_with_unmapped) {
		$pm->start and next;
		chr_map_stat($bam, $part, $samtools);	
		$pm->finish;	
	}
	$pm->wait_all_children;
	print "Finished all chrs mapping stat!\n";
	%bam_info = %{summary_map("$outDir/$sampleId.chr.stat.txt")};
	if ($target) {
		my $flank_shell=<<FLANK;
awk '{print \$1"\\t"\$2-$interval_padding"\\t"\$2-1"\\n"\$1"\\t"\$3+1"\\t"\$3+$interval_padding}' $target | sort -k 1,1 -k 2,2n  > $outDir/near.bed
$bedtools merge -i $outDir/near.bed > $outDir/near.merge.bed
$bedtools subtract -a $outDir/near.merge.bed -b $target > $outDir/flank.bed
rm $outDir/near*
FLANK
		system("$flank_shell") == 0 || die $!;
		my $chr_ref = region_split($target, "$outDir/split_region");
		region_split("$outDir/flank.bed", "$outDir/split_region");
		print "Start depth stat by chr!\n";
		system(">$outDir/$sampleId.chr.depth.txt") == 0 || die $!;
		foreach my $task (@$chr_ref) {
			foreach my $input ($target, "$outDir/flank.bed") {
				my $prefix = basename($input);
			        $prefix =~ s/.bed//;
				$pm->start and next;
				chr_depth($bam, "$outDir/split_region/$prefix.$task.bed", $samtools);
				$pm->finish;
			}
		}
		$pm->wait_all_children;
		print "Finished all chrs depth stat!\n";
		%depth_info = %{summary_depth("$outDir/$sampleId.chr.depth.txt")};
		system("rm -rf $outDir/flank.bed $outDir/split_region") == 0 || die $!;
		my $map_rate = sprintf("%.2f\%", $bam_info{'map'}/$bam_info{'reads'}*100);
		my $uniq_rate = sprintf("%.2f\%", $bam_info{'uniq'}/$bam_info{'reads'}*100); 
		my $dup_rate = sprintf("%.2f\%", $bam_info{'dup'}/$bam_info{'reads'}*100);
		my $mis_rate = sprintf("%.2f\%", $bam_info{'mismatch'}/$bam_info{'bases'}*100);
		my $capture_rate = sprintf("%.2f\%", $depth_info{'target'}{'reads'}/$bam_info{'map'}*100);
		my $target_sequence_rate = sprintf("%.2f\%", $depth_info{'target'}{'sequences'}/$bam_info{'effective_bases'}*100); 
		my $flank_sequence_rate = sprintf("%.2f\%", $depth_info{'flank'}{'sequences'}/$bam_info{'effective_bases'}*100); 
		my $ave_target_depth = sprintf("%.2f", $depth_info{'target'}{'sequences'}/$depth_info{'target'}{'bases'});
		my $ave_flank_depth = sprintf("%.2f", $depth_info{'flank'}{'sequences'}/$depth_info{'flank'}{'bases'});
		my $target_cov_rate = sprintf("%.2f\%", $depth_info{'target'}{'cov'}/$depth_info{'target'}{'bases'}*100); 
		my $flank_cov_rate = sprintf("%.2f\%", $depth_info{'flank'}{'cov'}/$depth_info{'flank'}{'bases'}*100); 
		my $target_4x_rate = sprintf("%.2f\%", $depth_info{'target'}{'x4'}/$depth_info{'target'}{'bases'}*100);
		my $target_10x_rate = sprintf("%.2f\%", $depth_info{'target'}{'x10'}/$depth_info{'target'}{'bases'}*100);
		my $target_20x_rate = sprintf("%.2f\%", $depth_info{'target'}{'x20'}/$depth_info{'target'}{'bases'}*100);
		my $target_30x_rate = sprintf("%.2f\%", $depth_info{'target'}{'x30'}/$depth_info{'target'}{'bases'}*100);
		my $flank_4x_rate = sprintf("%.2f\%", $depth_info{'flank'}{'x4'}/$depth_info{'flank'}{'bases'}*100);
		my $flank_10x_rate = sprintf("%.2f\%", $depth_info{'flank'}{'x10'}/$depth_info{'flank'}{'bases'}*100);
		my $flank_20x_rate = sprintf("%.2f\%", $depth_info{'flank'}{'x20'}/$depth_info{'flank'}{'bases'}*100);
		my $flank_30x_rate = sprintf("%.2f\%", $depth_info{'flank'}{'x30'}/$depth_info{'flank'}{'bases'}*100);
		$output=<<OUTPUT;
Sample\t$sampleId
Reads_number\t$bam_info{'reads'}
Bases_number\t$bam_info{'bases'}
Total_effective_bases\t$bam_info{'effective_bases'}
Reads_mapped_genome\t$bam_info{'map'}
Reads_mapped_genome_rate\t$map_rate
Reads_mapped_target\t$depth_info{'target'}{'reads'}
Reads_capture_rate\t$capture_rate
Reads_mapped_flank\t$depth_info{'flank'}{'reads'}
Uniq_mapped_rate\t$uniq_rate
Duplication_rate\t$dup_rate
Mismatch_bases_rate\t$mis_rate
Target_bases\t$depth_info{'target'}{'bases'}
Flank_bases\t$depth_info{'flank'}{'bases'}
Effective_bases_on_target\t$depth_info{'target'}{'sequences'}
Effective_bases_on_flank\t$depth_info{'flank'}{'sequences'}
Fraction_effective_bases_on_target\t$target_sequence_rate
Fraction_effective_bases_on_flank\t$flank_sequence_rate
Average_effective_depth_on_target\t$ave_target_depth
Average_effective_depth_on_flank\t$ave_flank_depth
Base_covered_on_target\t$depth_info{'target'}{'cov'}
Base_covered_on_flank\t$depth_info{'flank'}{'cov'}
Target_coverage_rate\t$target_cov_rate
Flank_coverage_rate\t$flank_cov_rate
Fraction_target_covered_at_least_4x\t$target_4x_rate
Fraction_target_covered_at_least_10x\t$target_10x_rate
Fraction_target_covered_at_least_20x\t$target_20x_rate
Fraction_target_covered_at_least_30x\t$target_30x_rate
Fraction_flank_covered_at_least_4x\t$flank_4x_rate
Fraction_flank_covered_at_least_10x\t$flank_10x_rate
Fraction_flank_covered_at_least_20x\t$flank_20x_rate
Fraction_flank_covered_at_least_30x\t$flank_30x_rate
OUTPUT
	} else {
		print "Start depth stat by chr!\n";
		system(">$outDir/$sampleId.chr.depth.txt") == 0 || die $!;
		foreach my $task (@all_chr) {
			$pm->start and next;
			chr_depth($bam, $task, $samtools);
			$pm->finish;
		}
		$pm->wait_all_children;
		print "Finished all chrs depth stat!\n";
		%depth_info = %{summary_depth("$outDir/$sampleId.chr.depth.txt")};
		my $map_rate = sprintf("%.2f\%", $bam_info{'map'}/$bam_info{'reads'}*100);
		my $uniq_rate = sprintf("%.2f\%", $bam_info{'uniq'}/$bam_info{'reads'}*100);
		my $dup_rate = sprintf("%.2f\%", $bam_info{'dup'}/$bam_info{'reads'}*100);
		my $mis_rate = sprintf("%.2f\%", $bam_info{'mismatch'}/$bam_info{'bases'}*100);
		my $genome_sequence_rate = sprintf("%.2f", $depth_info{$sampleId}{'sequences'}/$bam_info{'effective_bases'}*100); 
		my $ave_genome_depth = sprintf("%.2f", $depth_info{$sampleId}{'sequences'}/$depth_info{$sampleId}{'bases'});
		my $genome_cov_rate = sprintf("%.2f\%", $depth_info{$sampleId}{'cov'}/$depth_info{$sampleId}{'bases'}*100); 
		my $genome_4x_rate = sprintf("%.2f\%", $depth_info{$sampleId}{'x4'}/$depth_info{$sampleId}{'bases'}*100);
		my $genome_10x_rate = sprintf("%.2f\%", $depth_info{$sampleId}{'x10'}/$depth_info{$sampleId}{'bases'}*100);
		my $genome_20x_rate = sprintf("%.2f\%", $depth_info{$sampleId}{'x20'}/$depth_info{$sampleId}{'bases'}*100);
		my $genome_30x_rate = sprintf("%.2f\%", $depth_info{$sampleId}{'x30'}/$depth_info{$sampleId}{'bases'}*100);
		$output=<<OUTPUT;
Sample\t$sampleId
Reads_number\t$bam_info{'reads'}
Bases_number\t$bam_info{'bases'}
Total_effective_bases\t$bam_info{'effective_bases'}
Reads_mapped_genome\t$bam_info{'map'}
Reads_mapped_genome_rate\t$map_rate
Uniq_mapped_rate\t$uniq_rate
Mismatch_bases_rate\t$mis_rate
Duplication_rate\t$dup_rate
Genome_sizes\t$depth_info{$sampleId}{'bases'}
Average_effective_depth_on_genome\t$ave_genome_depth
Base_covered_on_genome\t$depth_info{$sampleId}{'cov'}
Genome_coverage_rate\t$genome_cov_rate
Fraction_genome_covered_at_least_4x\t$genome_4x_rate
Fraction_genome_covered_at_least_10x\t$genome_10x_rate
Fraction_genome_covered_at_least_20x\t$genome_20x_rate
Fraction_genome_covered_at_least_30x\t$genome_30x_rate
OUTPUT
	}
	open R, ">$outDir/$sampleId.bam.stat.txt" or die $!;
	print R $output;
	close R;
}

sub region_split {
	my $file = shift;
	my $outDir = shift;
	my %chr_region;
	my @chrs; 
	my $prefix = basename($file);
	$prefix =~ s/.bed//;
	open F, $file or die $!;
	while (<F>) {
		next if (/#/);
		next if (/^\s*$/);
		my @arr = split /\s+/;
		$chr_region{$arr[0]} .= $_;
	}
	close F;
	@chrs = sort {$a cmp $b} keys %chr_region;
	system("mkdir -p $outDir") == 0 || die $!;
	foreach my $chr (@chrs) {
		open OUT, ">$outDir/$prefix.$chr.bed" or die $!;
		print OUT $chr_region{$chr};
		close OUT;
	}
	return \@chrs;
}

sub chr_depth {
	my $bam = shift;
	my $target = shift;
	my $samtools = shift;
	my $outDir = dirname($bam);
	my $sampleId = basename($bam);
	$sampleId =~ s/.bam//;
	$sampleId =~ s/.final//;
	my $bases = 0;
	my $sequences = 0;
	my $cov = 0;
	my $x4 = 0;
	my $x10 = 0;
	my $x20 = 0;
	my $x30 = 0;
	my $reads = 0;
	my ($prefix, $chr);
	if (-e $target) {
		my @arr = split(/\./, basename($target));
		$prefix = $arr[-3];
		$prefix = 'target' if ($prefix ne 'flank');
		$chr = $arr[-2];
		open IN, "$samtools depth -a -d 0 -b $target $bam |" or die $!;
		open SR, "$samtools view -L $target -F 4 $bam |" or die $!;
	} else {
		$prefix = $sampleId;
		$chr = $target;
		open IN, "$samtools depth -a -d 0 -r $chr $bam |" or die $!;
		open SR, "$samtools view -F 4 $bam $chr |" or die $!;
	}
	while (<IN>) {
		my @arr = split /\t/;
		$bases++;
		$cov++ if ($arr[2] > 0);
		$x4++ if ($arr[2] > 3);
		$x10++ if ($arr[2] > 9);
		$x20++ if ($arr[2] > 19);
		$x30++ if ($arr[2] > 29);
		$sequences += $arr[2];
	}
	close IN;

	while (<SR>) {
		$reads++;	
	}
	close SR;
	open STAT, ">>$outDir/$sampleId.chr.depth.txt" or die $!;
	print STAT "$prefix\t$chr\t$bases\t$cov\t$x4\t$x10\t$x20\t$x30\t$sequences\t$reads\n";
	close STAT;
}

sub summary_depth {
	my $file = shift;
	my %info;
	open IN, $file or die $!;
	while (<IN>) {
		my @arr = split /\t/;
		$info{$arr[0]}{'bases'} += $arr[2];
		$info{$arr[0]}{'cov'} += $arr[3];
		$info{$arr[0]}{'x4'} += $arr[4];
		$info{$arr[0]}{'x10'} += $arr[5];
		$info{$arr[0]}{'x20'} += $arr[6];
		$info{$arr[0]}{'x30'} += $arr[7];
		$info{$arr[0]}{'sequences'} += $arr[8];
		$info{$arr[0]}{'reads'} += $arr[9];
	} 
	close IN;
	return \%info;
}

sub bam_stat {
	my $bam = shift;
	my $samtools =shift;
	my %info;
	open BAM, "$samtools view $bam |" or die $!;
	while (<BAM>) {
		my @arr = split /\t/;
		$info{'reads'}++;
		$info{'bases'} += length($arr[9]);
		$info{'map'}++ unless ($arr[1] & 4);
		$info{'dup'}++ if ($arr[1] & 1024);
	}
	close BAM;
	return \%info;
}

sub chr_map_stat {
	$bam = shift;
	$chr = shift;
	$samtools = shift;
	my $outDir = dirname($bam);
	my $sampleId = basename($bam);
	$sampleId =~ s/.bam//;
	$sampleId =~ s/.final//;
	my $reads = 0;
	my $bases = 0;
	my $map = 0;
	my $uniq_map_reads = 0;
	my $dup = 0;
	my $effective_bases = 0;
	my $mismatch_bases = 0;
	if ($chr ne 'unmapped') {
		open BAM, "$samtools view $bam $chr |" or die $!;
	} else {
		open BAM, "$samtools view -f 4 $bam |" or die $!;
	}
	while (<BAM>) {
		my $line = $_;
		my @arr = split /\t/;
		if ($chr eq 'unmapped') {
			next unless ($arr[2] eq '*');
		}
		unless ($arr[1] & 2048) {
			$reads++;
			$bases += length($arr[9]);
			$map++ unless ($arr[1] & 4);
			$uniq_map_reads++ if ($arr[4] > 1);
			$dup++ if ($arr[1] & 1024);
			unless ($arr[1] & 4) {
				unless ($arr[1] & 1024) {
					$effective_bases += length($arr[9]);
				}
			}
			unless ($arr[1] & 4) {
				$line =~ m/\s+(MD\S+)\s+/;
				my @md = split ":", $1;
				my $str = $md[-1];
				$str =~ s/\^[A|T|C|G]//gi;
				my $tmp_count = $str =~ tr/ATCGatcg//;
				$mismatch_bases += $tmp_count;
			}
		}
	}
	close BAM;	
	open OUT, ">>$outDir/$sampleId.chr.stat.txt" or die $!;
	print OUT "$chr\t$reads\t$bases\t$map\t$uniq_map_reads\t$dup\t$effective_bases\t$mismatch_bases\n";
	close OUT;
}

sub summary_map {
	my $file = shift;
	my %info;
	open IN, $file or die $!;
	while (<IN>) {
		my @arr = split /\t/;
		$info{'reads'} += $arr[1];
		$info{'bases'} += $arr[2];
		$info{'map'} += $arr[3];
		$info{'uniq'} += $arr[4];
		$info{'dup'} += $arr[5];
		$info{'effective_bases'} += $arr[6];
		$info{'mismatch'} += $arr[7];
	} 
	close IN;
	return \%info;
}

1;
