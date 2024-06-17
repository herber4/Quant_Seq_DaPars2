# Differential APA analysis from Quant Seq Rev data using DaPars2

This is an imperfect protocol and still under development

## Align reads with STAR, soft alignment constraints. Convert sam to bam, sort, index and gather coverage info

`
ml star/2.7.10a

for i in *_R1_001.fastq.gz; do name=$(basename ${i} _R1_001.fastq.gz);

STAR --runThreadN 32 \
--readFilesCommand zcat --genomeDir /data2/lackey_lab/austin/dbs/star \
--readFilesIn ${name}_R1_001.fastq.gz ${name}_R2_001.fastq.gz \
--outFilterType BySJout --outFilterMultimapNmax 200 \
--outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 \
--alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.6 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
--limitOutSJcollapsed 5000000 \
--quantMode TranscriptomeSAM \
--outSAMattributes NH HI NM MD --outSAMtype BAM Unsorted \
--limitBAMsortRAM 2000000000 \
--outFileNamePrefix ${name}_ ;
done

module load samtools/1.10
for B in *_Aligned.out.bam; do
    N=$(basename $B _Aligned.out.bam) ;
    samtools view --threads 8 -bS $B | samtools sort -O BAM -o ${N}.sorted.bam ;
done

for b in *.sorted.bam ; do
	samtools index $b ;
done

for b in *.sorted.bam; do
    N=$(basename $B .sorted.bam) ;
    samtools coverage -o ${N}.txt $b ;
done
`

## Generate genome coverage bedgraphs from bam files

`
for b in *.sorted.bam; do
	N=$(basename $b .sorted.bam) ;
	bedtools genomecov -ibam $b -bga -split -trackline > $N.bedgraph ;
done
`

## Generate samtools coverages to coverage input in DaPars2

`
module load samtools/1.10

for b in *.sorted.bam; do
    N=$(basename $b .sorted.bam) ;
    samtools coverage -o ${N}.txt $b ;
done
`
# run the reduce_chrs.sh script to get chromosomes from samtools coverage output

`
#!/bin/bash

for t in *.txt ; do
	N=$(basename $t .txt) ;
	head -n 26 $t > $N.reduced.txt ;
done
`

# run this script to get coverage sums from each sample

`
#!/bin/bash

for t in *.txt ; do
	N=$(basename $t .txt) ;
	awk 'NR > 1 { sum += $4 } END { print sum }' $t > $N.counts ;
done

cat *.counts > merged.counts
`

## Running Dapars2

- examples of files you need are in the github_dapars2 dir
- make sure all bedgrpahs and necessary files are in the same directory to run
- items needed to run:
bedgraph for each sample
chrList.txt
dapars.config
sample_coverages.txt
RefSeq_hg38_3UTR_annotation.bed

`
conda activate dapars

python3 ../DaPars2/src/DaPars2_Multi_Sample_Multi_Chr.py dapars.config chrList.txt
`

# Moving raw dapars2 for statistical analysis in R

You need to move all chromosome output into one dir, as well as create a meta data file like the sample_meta.txt, you can do this in excel

`
mkdir output
cp */*.txt output/
`
