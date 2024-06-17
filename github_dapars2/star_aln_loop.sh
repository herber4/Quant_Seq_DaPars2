#!/bin/bash
#SBATCH --job-name=star
#SBATCH -n 32
#SBATCH --partition=compute
#SBATCH --time=100:00:00
#SBATCH --mem=128G
#SBATCH --mail-type=all
#SBATCH --output=/data2/lackey_lab/DownloadedSequenceData/austin/dapars/star.%j.out
#SBATCH --error=/data2/lackey_lab/DownloadedSequenceData/austin/dapars/star.%j.err
#SBATCH --mail-user=herber4@clemson.edu
# Load change dir
cd /data2/lackey_lab/DownloadedSequenceData/austin/dapars/fqs

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

for b in *.sorted.bam; do
    N=$(basename $B .sorted.bam) ;
    samtools coverage -o ${N}.txt $b ;
done
