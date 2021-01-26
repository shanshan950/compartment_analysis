# compartment_analysis testing
Here is an example to show the 250kb resolution compartment analysis.
You will need fragment cis reads file: <frag_1> <frag_2> <reads>
hg19.genome_split_250k file format is <chr1> <start> <end> <bin1> 
hg19.simplify_chr.fra.bin.250k <frag_1> <chr1_bin1>
H3K4me3 peaks bed file: <chr> <start> <end> <peak_name>
```
# initialize variables
ref=ref/
lib=lib/
frag_loop=
output=./
bin=250k # resolution e.g. 250k; 500k; 1M
chipseq_peak_bedfile=$4

# create directory for output
mkdir $output
mkdir $output/$bin

# the output is for identifying the (-/+) of PC value
bedtools intersect -wa -wb -a $ref/hg19.genome_split_250k -b $chipseq_peak_bedfile > $output/$bin/chip.file

# convert frag_to_frag files to bin_to_bin files, $frag_loop should be sorted
$lib/1.frag.to.bin.merge.from.loop_V2.pl $ref/hg19.simplify_chr.fra.bin.$bin $frag_loop > $output/$bin/$bin.sparse.matrix

# create matrix for each chromosome matrix
for i in {1..22} X Y;do
        cat $output/$bin/$bin.sparse.matrix | grep chr${i}_bin >> $output/$bin/$bin.sparse.matrix/chr${i}.matrix
done

cd $output/$bin

# distance correction and calculate PCA
$lib/3.generate.chr.matrix_component.r $output $bin

# summarize H3K4me3 peak counts to define negative PC1 bins and positive PC1 bins, 
for i in {1..22} X Y;do
        $lib/compartment_calling.pl chip.file chr${i}.matrix.component > chr${i}.bin.$bin.chip.freq.score
done
# plot heatmaps with the first 3 PC values, you may need to check if PC1 is the best for explaining compartment for each chromosome

$lib/5.plot.three.component.r `pwd` $bin $ref

# if thre is any chromosome explained better by PC2, you need to change the default "1" to "2" for the chromosome in $lib/6.plot.one.component.r
$lib/6.plot.one.component.r $output $bin `pwd` $ref

# the output contains <bin> <compartmrnt A/B> <PC values>
cat chr*bin.$bin.chip.freq.score.A.B.label | grep -v "chrX" | grep -v "chrY" > $output.$bin.merge.label
```
