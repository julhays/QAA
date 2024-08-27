#!/bin/bash

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH -c 1                 #optional: number of cpus, default is 1
#SBATCH --mem=32GB                        #optional: amount of memory, default is 4GB per cpu
#SBATCH --time=0-8
#SBATCH --mail-user=jkhay@uoregon.edu     #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --job-name=cutadapt           #optional: job name
#SBATCH --output=cutadapt_%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=cutadapt_%j.err        #optional: file to store stderr from job, %j adds the assigned jobID


conda activate QAA

/usr/bin/time -v cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o cutadapt_out/16_R1_trimmed.fastq \
    -p cutadapt_out/16_R2_trimmed.fastq \
    /projects/bgmp/shared/2017_sequencing/demultiplexed/16_3D_mbnl_S12_L008_R1_001.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/16_3D_mbnl_S12_L008_R2_001.fastq.gz
/usr/bin/time -v cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o cutadapt_out/6_R1_trimmed.fastq \
    -p cutadapt_out/6_R2_trimmed.fastq \
    /projects/bgmp/shared/2017_sequencing/demultiplexed/6_2D_mbnl_S5_L008_R1_001.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/6_2D_mbnl_S5_L008_R2_001.fastq.gz

exit





