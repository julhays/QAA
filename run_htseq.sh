#!/bin/bash

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH -c 8                 #optional: number of cpus, default is 1
#SBATCH --mem=32GB                        #optional: amount of memory, default is 4GB per cpu
#SBATCH --time=0-8
#SBATCH --mail-user=jkhay@uoregon.edu     #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --job-name=htseq           #optional: job name
#SBATCH --output=htseq_%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=htseq_%j.err        #optional: file to store stderr from job, %j adds the assigned jobID


conda activate QAA

/usr/bin/time -v htseq-count -f sam -r name --stranded=yes -c 6_stranded.genecount Mus_aligned_out/6_Mus_Aligned.out.sam mmus_files/Mus_musculus.GRCm39.112.gtf
/usr/bin/time -v htseq-count -f sam -r name --stranded=yes -c 16_stranded.genecount Mus_aligned_out/16_Mus_Aligned.out.sam mmus_files/Mus_musculus.GRCm39.112.gtf
/usr/bin/time -v htseq-count -f sam -r name --stranded=reverse -c 6_reverse.genecount Mus_aligned_out/6_Mus_Aligned.out.sam mmus_files/Mus_musculus.GRCm39.112.gtf
/usr/bin/time -v htseq-count -f sam -r name --stranded=reverse -c 16_reverse.genecount Mus_aligned_out/16_Mus_Aligned.out.sam mmus_files/Mus_musculus.GRCm39.112.gtf

exit



