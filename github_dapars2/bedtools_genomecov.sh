#!/bin/bash
#SBATCH --job-name=bedtools
#SBATCH -n 32
#SBATCH --partition=compute
#SBATCH --time=100:00:00
#SBATCH --mem=128G
#SBATCH --mail-type=all
#SBATCH --output=/data2/lackey_lab/DownloadedSequenceData/austin/dapars/genomecov.%j.out
#SBATCH --error=/data2/lackey_lab/DownloadedSequenceData/austin/dapars/genomecov.%j.err
#SBATCH --mail-user=herber4@clemson.edu
# Load change dir
cd /data2/lackey_lab/DownloadedSequenceData/austin/dapars/fqs
conda activate dapars

for b in *.sorted.bam; do
	N=$(basename $b .sorted.bam) ;
	bedtools genomecov -ibam $b -bga -split -trackline > $N.bedgraph ;
done
