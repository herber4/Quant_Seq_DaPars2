#!/bin/bash
#SBATCH --job-name=star
#SBATCH -n 16
#SBATCH --partition=compute
#SBATCH --time=100:00:00
#SBATCH --mem=64G
#SBATCH --mail-type=all
#SBATCH --output=/data2/lackey_lab/DownloadedSequenceData/austin/dapars/sam_covs.%j.out
#SBATCH --error=/data2/lackey_lab/DownloadedSequenceData/austin/dapars/sam_covs.%j.err
#SBATCH --mail-user=herber4@clemson.edu
# Load change dir
cd /data2/lackey_lab/DownloadedSequenceData/austin/dapars/fqs

module load samtools/1.10

for b in *.sorted.bam; do
    N=$(basename $b .sorted.bam) ;
    samtools coverage -o ${N}.txt $b ;
done
