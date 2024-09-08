#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH -c 8
#SBATCH --mem=100G
#SBATCH --time=0-3
#SBATCH --job-name STAR_align
#SBATCH --output=staralign_%j.out
#SBATCH --output=staralign_%j.err
#SBATCH --mail-user=jkhay@uoregon.edu
#SBATCH --mail-type=ALL

conda activate QAA

/usr/bin/time -v STAR --runThreadN 8 --runMode alignReads \
    --outFilterMultimapNmax 3 \
    --outSAMunmapped Within KeepPairs \
    --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
    --readFilesCommand zcat \
    --readFilesIn trim_out/6_R1_paired_filtered.fastq.gz trim_out/6_R2_paired_filtered.fastq.gz \
    --genomeDir Mus_musculus.GRCm39.dna.ens112.STAR_2.7.11b \
    --outFileNamePrefix 6_Mus_

/usr/bin/time -v STAR --runThreadN 8 --runMode alignReads \
    --outFilterMultimapNmax 3 \
    --outSAMunmapped Within KeepPairs \
    --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
    --readFilesCommand zcat \
    --readFilesIn trim_out/16_R1_paired_filtered.fastq.gz trim_out/16_R2_paired_filtered.fastq.gz \
    --genomeDir Mus_musculus.GRCm39.dna.ens112.STAR_2.7.11b \
    --outFileNamePrefix 16_Mus_

exit