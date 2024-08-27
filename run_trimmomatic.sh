#!/bin/bash

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH -c 8                 #optional: number of cpus, default is 1
#SBATCH --mem=32GB                        #optional: amount of memory, default is 4GB per cpu
#SBATCH --time=0-8
#SBATCH --mail-user=jkhay@uoregon.edu     #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --job-name=trimmomatic           #optional: job name
#SBATCH --output=trimmomatic_%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=trimmomatic_%j.err        #optional: file to store stderr from job, %j adds the assigned jobID


conda activate QAA

/usr/bin/time -v trimmomatic PE \
    cutadapt_out/16_R1_trimmed.fastq cutadapt_out/16_R2_trimmed.fastq \
    trim_out/16_R1_paired_filtered.fastq.gz \
    trim_out/16_R1_unpaired_filtered.fastq.gz \
    trim_out/16_R2_paired_filtered.fastq.gz \
    trim_out/16_R2_unpaired_filtered.fastq.gz \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:5:15 \
    MINLEN:35

/usr/bin/time -v trimmomatic PE \
    cutadapt_out/6_R1_trimmed.fastq cutadapt_out/6_R2_trimmed.fastq \
    trim_out/6_R1_paired_filtered.fastq.gz \
    trim_out/6_R1_unpaired_filtered.fastq.gz \
    trim_out/6_R2_paired_filtered.fastq.gz \
    trim_out/6_R2_unpaired_filtered.fastq.gz \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:5:15 \
    MINLEN:35

exit



