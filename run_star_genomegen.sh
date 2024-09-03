#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH -c 8
#SBATCH --mem=100G
#SBATCH --time=0-3
#SBATCH --job-name STAR_index
#SBATCH --output=star_%j.out
#SBATCH --output=star_%j.err
#SBATCH --mail-user=jkhay@uoregon.edu
#SBATCH --mail-type=ALL

conda activate QAA

/usr/bin/time -v STAR \
    --runThreadN 8 --runMode genomeGenerate \
    --genomeDir Mus_musculus.GRCm39.dna.ens112.STAR_2.7.11b \
    --genomeFastaFiles mmus_files/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa \
    --sjdbGTFfile mmus_files/Mus_musculus.GRCm39.112.gtf

exit