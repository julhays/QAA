#!/bin/bash

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH -c 8                 #optional: number of cpus, default is 1
#SBATCH --mem=100GB                        #optional: amount of memory, default is 4GB per cpu
#SBATCH --time=0-8
#SBATCH --mail-user=jkhay@uoregon.edu     #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --job-name=readlen_dist           #optional: job name
#SBATCH --output=readlen_dist_%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=readlen_dist_%j.err        #optional: file to store stderr from job, %j adds the assigned jobID


conda activate QAA

/usr/bin/time -v ./read_len_dist.py \
    -p1 'trim_out/16_R1_paired_filtered.fastq.gz' \
    -u1 'trim_out/16_R1_unpaired_filtered.fastq.gz' \
    -p2 'trim_out/16_R2_paired_filtered.fastq.gz' \
    -u2 'trim_out/16_R2_unpaired_filtered.fastq.gz' \
    -o 16

/usr/bin/time -v ./read_len_dist.py \
    -p1 'trim_out/6_R1_paired_filtered.fastq.gz' \
    -u1 'trim_out/6_R1_unpaired_filtered.fastq.gz' \
    -p2 'trim_out/6_R2_paired_filtered.fastq.gz' \
    -u2 'trim_out/6_R2_unpaired_filtered.fastq.gz' \
    -o 6

exit





