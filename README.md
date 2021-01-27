# compartment_analysis testing
Here is an example to show the 250kb resolution compartment analysis.
```
# data requirement:
# fragment cis reads file: <frag_1> <frag_2> <reads>
# hg19.genome_split_250k file format is <chr> <start> <end> <bin> 
# hg19.simplify_chr.fra.bin.250k <frag> <chr_bin>
# H3K4me3 peaks bed file: <chr> <start> <end> <peak_name>

# initialize variables
ref=ref/
lib=lib/
frag_loop= # example data 
output=test
bin=250k # resolution e.g. 250k; 500k; 1M
binbed=ref/hg19.genome_split_250k
chipsePeakbed=$4 # H3K4me3 peaks bed file

# create directory for output
mkdir $output
mkdir $output/$bin

# convert frag_to_frag files to bin_to_bin files, $frag_loop should be sorted
# eHiC fragment with "+" or "-" following frag ID
cat $frag_loop | sed s/"+"//g | sed s/"-"//g | $lib/1.frag.to.bin.merge.from.loop_V2.pl $ref/hg19.simplify_chr.fra.bin.$bin - > $output/$bin/$bin.sparse.matrix
# conventional HiC with no "+" or "-"
cat $frag_loop | $lib/1.frag.to.bin.merge.from.loop_V2.pl $ref/hg19.simplify_chr.fra.bin.$bin - > $output/$bin/$bin.sparse.matrix

# create matrix for each chromosome matrix
cd $output/$bin/
cat $bin.sparse.matrix | sed s/"_"/\\t/g | awk '{print $1"_"$2 "\t" $3"_"$4 "\t" $5 >$1".matrix"}'

# distance correction and calculate PCA
$lib/3.generate.chr.matrix_component.r ./ $binbed

# summarize H3K4me3 peak counts to define negative PC1 bins and positive PC1 bins
# the output is for identifying the (-/+) of PC value
for i in {1..22} X Y;do
        bedtools intersect -wa -wb -a $binbed -b $chipsePeakbed | $lib/compartment_calling.pl - chr${i}.matrix.component > chr${i}.bin.$bin.chip.freq.score
done

# plot heatmaps with the first 3 PC values, you may need to check if PC1 is the best for explaining compartment for each chromosome
$lib/5.plot.three.component.r ./ $binbed 250k

# Visually check 3 PC compartments with each heatmap. In most cases, PC1 is the best matched. 

```
