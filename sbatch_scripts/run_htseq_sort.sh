#!/bin/bash

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH -c 8                 #optional: number of cpus, default is 1
#SBATCH --mem=32GB                        #optional: amount of memory, default is 4GB per cpu
#SBATCH --time=0-8
#SBATCH --mail-user=jkhay@uoregon.edu     #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --job-name=htseq_sort           #optional: job name
#SBATCH --output=htseqsort_%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=htseqsort_%j.err        #optional: file to store stderr from job, %j adds the assigned jobID


conda activate QAA

/usr/bin/time -v htseq-count -r name --stranded=yes Mus_aligned_out/sorted_6_Mus_Aligned.out.sam mmus_files/Mus_musculus.GRCm39.112.gtf > counts_files/sort_6_stranded.genecount
/usr/bin/time -v htseq-count -r name --stranded=yes Mus_aligned_out/sorted_16_Mus_Aligned.out.sam mmus_files/Mus_musculus.GRCm39.112.gtf > counts_files/sort_16_stranded.genecount
/usr/bin/time -v htseq-count -r name --stranded=reverse Mus_aligned_out/sorted_6_Mus_Aligned.out.sam mmus_files/Mus_musculus.GRCm39.112.gtf > counts_files/sort_6_reverse.genecount
/usr/bin/time -v htseq-count -r name --stranded=reverse Mus_aligned_out/sorted_16_Mus_Aligned.out.sam mmus_files/Mus_musculus.GRCm39.112.gtf > counts_files/sort_16_reverse.genecount

exit



