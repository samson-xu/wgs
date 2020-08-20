package CNV;

use File::Basename;
use WriteShell;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(cnvkit);

sub cnvkit {
	my $cnv_soft = 'cnvkit.py';
	my $outDir = shift;
	my $target = shift;
	my $thread = shift;
	my $bam = shift;
	my $ref = shift;
	my $access = shift;
	my $filter = shift;
	my $annot = shift;
	my $flag = shift;
	my $cnv_db = shift;
	my $method;
	if ($target) {
		$method = 'hybrid';
		$target = "--targets $target";
	} else {
		$method = 'wgs';
		$target = "";
	}
	my $cnv_shell = "# pipeline for detect cnv\n";
	unless($cnv_db) {
		my $dbDir = "$outDir/db";
		my $db_shell=<<DB;
# shell for create cnv db 
# 初始化conda环境
source activate
conda deactivate

# 切换cnvkit环境
conda activate cnvkit-0.9.6
$cnv_soft batch \\
    --method $method \\
    --processes $thread \\
    --normal $bam \\
    --fasta $ref \\
    --access $access \\
    --output-reference $dbDir/ref.cnn \\
    --output-dir $dbDir $target

#$cnv_soft reference $dbDir/*.cnn --fasta $ref --output $dbDir/ref.cnn

$cnv_soft sex $dbDir/*.targetcoverage.cnn > $dbDir/sample.sex.info
conda deactivate

DB
	write_shell($db_shell, "$dbDir/build_ref.sh");	
	$cnv_db = "$dbDir/ref.cnn";
	$cnv_shell .= "$dbDir/build_ref.sh >$dbDir/build_ref.sh.o 2>$dbDir/build_ref.sh.e\n";
	}
	my $callDir = "$outDir/cnv";
	my $call_shell=<<CALL;
# shell for detect cnv
# 初始化conda环境      
source activate        
conda deactivate  

# 切换cnvkit环境       
conda activate cnvkit-0.9.6

$cnv_soft batch $bam \\
--method $method \\
--processes $thread \\
--reference $cnv_db \\
--output-dir $callDir

$cnv_soft sex $callDir/*.cnr > $callDir/sample.sex.txt

for i in \$(ls $callDir/*.cnr)
do {
        dir=`dirname \$i`
        file=`basename \$i`
        prefix=\${file%.cnr}
        $cnv_soft segmetrics --ci -s $callDir/\$prefix.cns -o $callDir/\$prefix.cns.seg \$i
        $cnv_soft call -t=-1.1,-0.4,0.3,0.7 $callDir/\$prefix.cns.seg --ploidy 2 -o $callDir/\$prefix.cns.call
        $filter $callDir/sample.sex.txt $callDir/\$prefix.cns.call > $callDir/\$prefix.cnv.bed
        awk '{print \$0\"\\t\"\$3-\$2}' $callDir/\$prefix.cnv.bed > $callDir/\$prefix.cnv.size
        mkdir -p \$dir/\$prefix
        $annot -SVinputFile $callDir/\$prefix.cnv.bed -svtBEDcol 5 -outputFile $callDir/\$prefix/\$prefix.cnv.annot
        mv \$dir/\$prefix/*.tsv $callDir
        rm -rf \$dir/\$prefix
        sed -i '1s/SV length\\s*SV type/SV length\\tCopy number\\tSV type/' $callDir/\$prefix.cnv.annot.tsv
}&
done
wait

$cnv_soft metrics $callDir/*.cnr -s $callDir/*.cns > $callDir/sample.metrics.txt
conda deactivate

CALL
	write_shell($call_shell, "$callDir/cnvkit.sh");	
	$cnv_shell .= "$callDir/cnvkit.sh >$callDir/cnvkit.sh.o 2>$callDir/cnvkit.sh.e\n";
	write_shell($cnv_shell, "$outDir/cnv.$flag.sh");
	
}
