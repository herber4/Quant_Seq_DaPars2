#!/bin/bash
#SBATCH --job-name=index
#SBATCH -n 8
#SBATCH --partition=compute
#SBATCH --time=100:00:00
#SBATCH --mem=32G
#SBATCH --mail-type=all
#SBATCH --output=/data2/lackey_lab/DownloadedSequenceData/austin/dapars/index.%j.out
#SBATCH --error=/data2/lackey_lab/DownloadedSequenceData/austin/dapars/index.%j.err
#SBATCH --mail-user=herber4@clemson.edu
# Load change dir
cd /data2/lackey_lab/DownloadedSequenceData/austin/dapars/fqs

module load samtools/1.10
for b in *.sorted.bam ; do
	samtools index $b ;
done
