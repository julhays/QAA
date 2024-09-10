#!/bin/bash

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH -c 1                 #optional: number of cpus, default is 1
#SBATCH --mem=32GB                        #optional: amount of memory, default is 4GB per cpu
#SBATCH --time=0-8
#SBATCH --mail-user=jkhay@uoregon.edu     #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --job-name=qual_dist           #optional: job name
#SBATCH --output=qual_dist_%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=qual_dist_%j.err        #optional: file to store stderr from job, %j adds the assigned jobID


conda activate bgmp_py312

/usr/bin/time -v ./qual_dist.py -f /projects/bgmp/shared/2017_sequencing/demultiplexed/16_3D_mbnl_S12_L008_R1_001.fastq.gz -n 101 -o '16_R1'
/usr/bin/time -v ./qual_dist.py -f /projects/bgmp/shared/2017_sequencing/demultiplexed/16_3D_mbnl_S12_L008_R2_001.fastq.gz -n 101 -o '16_R2'
/usr/bin/time -v ./qual_dist.py -f /projects/bgmp/shared/2017_sequencing/demultiplexed/6_2D_mbnl_S5_L008_R1_001.fastq.gz -n 101 -o '6_R1'
/usr/bin/time -v ./qual_dist.py -f /projects/bgmp/shared/2017_sequencing/demultiplexed/6_2D_mbnl_S5_L008_R2_001.fastq.gz -n 101 -o '6_R2'

exit





