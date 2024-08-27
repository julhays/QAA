# BI623 RNA-seq Quality Assessment Assignment - Lab Notebook
Jules Hays
## Due: 9/13/24

Python version: 3.12.5

Environments used: QAA (fastqc 0.12.1, cutadapt 4.9, trimmomatic 0.39, matplotlib 3.9.2)

__File directory__

Python scripts:
* Qual_dist.py - plots the quality score per base of reads (Part 1)
* read_len_dist.py - plots the read length dist of trimmed and filtered reads (Part 2)

sbatch scripts:
* 

---
### 8/24/24
Despite Leslie telling me not to be scared, I'm scared.

I logged into Talapas, navigated to my Bi623 directory, and cloned the repository. The repository is located at the following path:
```
/projects/bgmp/jkhay/bioinfo/Bi623/QAA
```
I then copied the Github QAA repo using my GitHub token. The repo is found at the following link: https://github.com/julhays/QAA

I started an interactive session on Talapas:
```
$ srun --account=bgmp --partition=bgmp --time=3:00:00 --pty bash
```

Next I am going to locate my data. The assignments for the demultiplexed file pair I will be working with is located in a file that I viewed with the following command:
```
cat /projects/bgmp/shared/Bi623/QAA_data_assignments.txt
```

Copying and pasting from this file, my assignments are the following:
```
Jules   16_3D_mbnl_S12_L008     6_2D_mbnl_S5_L008
```
These 2 files should each have a read 1 and read 2 of the index - a total of 4 files to analyze.

The files are located at the followng paths, with an arrow pointing to what I will refer to each graph as in my file naming:
```
/projects/bgmp/shared/2017_sequencing/demultiplexed/16_3D_mbnl_S12_L008_R1_001.fastq.gz -> 16_R1
/projects/bgmp/shared/2017_sequencing/demultiplexed/16_3D_mbnl_S12_L008_R2_001.fastq.gz -> 16_R2
/projects/bgmp/shared/2017_sequencing/demultiplexed/6_2D_mbnl_S5_L008_R1_001.fastq.gz -> 6_R1
/projects/bgmp/shared/2017_sequencing/demultiplexed/6_2D_mbnl_S5_L008_R2_001.fastq.gz -> 6_R2
```
As per the instruvtions, I will not move, copy or unzip the files and instead will reference them by path and use zcat. Lets perform some initial data exploration.

File size:
```
$ ls -lah /projects/bgmp/shared/2017_sequencing/demultiplexed/
-rw-r-----+ 1 coonrod is.racs.pirg.bgmp 468M Sep  1  2020 16_3D_mbnl_S12_L008_R1_001.fastq.gz
-rw-r-----+ 1 coonrod is.racs.pirg.bgmp 452M Aug 23  2017 16_3D_mbnl_S12_L008_R2_001.fastq.gz
-rw-r-----+ 1 coonrod is.racs.pirg.bgmp 571M Aug 23  2017 6_2D_mbnl_S5_L008_R1_001.fastq.gz
-rw-r-----+ 1 coonrod is.racs.pirg.bgmp 653M Aug 23  2017 6_2D_mbnl_S5_L008_R2_001.fastq.gz
```
Ok they are each about 400-600 Mb zipped.

Find number of lines:
```
$ zcat /projects/bgmp/shared/2017_sequencing/demultiplexed/16_3D_mbnl_S12_L008_R1_001.fastq.gz | wc -l
    32940788
$ zcat /projects/bgmp/shared/2017_sequencing/demultiplexed/16_3D_mbnl_S12_L008_R2_001.fastq.gz | wc -l
    32940788
$ zcat /projects/bgmp/shared/2017_sequencing/demultiplexed/6_2D_mbnl_S5_L008_R1_001.fastq.gz | wc -l
    44112976 
$ zcat /projects/bgmp/shared/2017_sequencing/demultiplexed/6_2D_mbnl_S5_L008_R2_001.fastq.gz | wc -l
    44112976
```
Divide these by 4 to get number of records:
```
16_3D_mbnl_S12_L008 - 8,235,197
6_2D_mbnl_S5_L008 - 11,028,244
```
The R1 and R2 files for each pair have the same number of lines/records, which is a good sign.

Read length:
```
$ zcat <filename> | sed -n '2~4p' | awk '{print length($0)}' | uniq
```
They all have a read length of 101 which is what was expected.

I also checked the Phred encoding on each file by running the following command for each file
```
zcat <filepath> | sed -n '4~4p' | head
```
Since all files contained symbols like '#' and '<', it is Phred+33 encoding.

Below is a summary table of my initial data analysis:

| File | size (M) | # lines | # records | read_len | Phred encoding |
|---|---|---|---|---|---|
| 16_R1 | 468 | 32,940,788 | 8,235,197 | 101 | Phred+33 |
| 16_R2 | 452 | 32,940,788 | 8,235,197 | 101 | Phred+33 |
| 6_R1 | 571 | 44,112,976  | 11,028,244 | 101 | Phred+33 |
| 6_R2 | 653 | 44,112,976  | 11,028,244 | 101 | Phred+33 |

The files look good to continue with quality analysis.


### Part 1 - Read quality score ditributions
Goal: Install and use FastQC for quality assessment of demultiplexed RNA-seq files, compare the quality assessments to those from my own software.

I started an interactive session on Talapas:
```
$ srun --account=bgmp --partition=bgmp --time=3:00:00 --pty bash
```

I made a new conda environment called ```QAA``` and installed FastQC into it with the following commands:
```
$ conda create --name QAA python=3.12
$ conda activate QAA
$ conda install fastqc
$ fastqc --version                    
FastQC v0.12.1
```

Time to use FastQC. I will use the following FastQc manual and tutorial for help:

https://home.cc.umanitoba.ca/~psgendb/doc/fastqc.help

https://hbctraining.github.io/Intro-to-ChIPseq/lessons/02_QC_FASTQC.html

First we will use FastQC to make graphs about each fastq file such as a per-base quality score distribution for R1 and R2 reads of each pair, and per-base N content for each file.

I will run fastqc on all 4 files:
```
/usr/bin/time -v fastqc -o fastqc_out -t 4 /projects/bgmp/shared/2017_sequencing/demultiplexed/16_3D_mbnl_S12_L008_R1_001.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/16_3D_mbnl_S12_L008_R2_001.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/6_2D_mbnl_S5_L008_R1_001.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/6_2D_mbnl_S5_L008_R2_001.fastq.gz
```
The -t specifies the number of cores which my source said to do 1 per file, which is why I did 4. -o specifies the name of the output directory.

It ran succesfully and took 2 minute 49 seconds and 99% CPU to run all 4 files. My fastqc outputs for each file are located in ```fastqc_out``` directory. I unziped each of the zip files into their own directories, 1 for each file named after the file, using the following commands:
```
$ unzip fastqc_out/16_3D_mbnl_S12_L008_R1_001_fastqc.zip
$ unzip fastqc_out/16_3D_mbnl_S12_L008_R2_001_fastqc.zip
$ unzip fastqc_out/6_2D_mbnl_S5_L008_R1_001_fastqc.zip 
$ unzip fastqc_out/6_2D_mbnl_S5_L008_R2_001_fastqc.zip 
```

Below I will paste the resulting plots for each file.

| File name | per_base_qual | per_base_n |
|---|---|---|
| 16_R1 | ![alt text](fastqc_plots/16_R1_per_base_qual.png) | ![alt text](fastqc_plots/16_R1_per_base_n.png) |
| 16_R2 | ![alt text](fastqc_plots/16_R2_per_base_qual.png) | ![alt text](fastqc_plots/16_R2_per_base_n.png) |
| 6_R1 | ![alt text](fastqc_plots/6_R1_per_base_qual.png) | ![alt text](fastqc_plots/6_R1_per_base_n.png) |
| 6_R2 | ![alt text](fastqc_plots/6_R2_per_base_qual.png) | ![alt text](fastqc_plots/6_R2_per_base_n.png) |

The per base quality score distribution is consistant with the per base n content plot in all 4 files because there's a small spike in n content at the beginning of the reads and it looks like the quality score of the beginning the reads is lower than the rest of the read.

Now, I will run these files through my own script ```qual_dist.py``` to calculate quality score distribution to see if it is similar to FastQC's distribution. To do this, I copied the python script I used to do this, plus the slurm script to run it, in the Demultiplex assignment using the below commands. I also needed bioinfo.py since my script used it for the ```convert_phred``` function.
```
$ cp /projects/bgmp/jkhay/bioinfo/Bi622/
Demultiplex/Assignment-the-first/qual_dist.py .
$ cp /projects/bgmp/jkhay/bioinfo/Bi622/Demultiplex/Assignment-the-first/R1_dist.sh .
$ cp /projects/bgmp/jkhay/bioinfo/Bi622/Demultiplex/Assignment-the-first/bioinfo.py .
```
I am going to read through these scripts to make sure no changes need to be made to suit my purposes. I changed the name of the R1_dist.sh to ```run_qual_dist.sh``` because I am going to use to to run all my files through at once since they aren't as big as the files we had in Demux. I changed one thing - my naming scheme for my graph output was different in the Demux assignment and I had just used a regex to get if the file was R1, R2, R3, or R4. Since this is no longer applicable i am just going to add another argparse parameter for the output file name and title of the distribution plot.

Time to run  ```run_qual_dist.sh```.
```
$ sbatch run_qual_dist.sh 
Submitted batch job 15875460
```
It succesffully ran, but it took a long time. If I do this again I will use 4 separate slurm scripts so I can run all 4 files in parallel. It took about 4-6 minutes per file for a total of 20 minutes and 8 seconds. Each run took about 99% CPU. The err and out files are in my ```Logs``` directory.

Below is a table to compare the per base average quality score of the 4 files:
| File name | fastq_per_base_qual | my_per_base_qual |
|---|---|---|
| 16_R1 | ![alt text](fastqc_plots/16_R1_per_base_qual.png) | ![alt text](qual_plots/16_R1_dist.png) |
| 16_R2 | ![alt text](fastqc_plots/16_R2_per_base_qual.png) | ![alt text](qual_plots/16_R2_dist.png) |
| 6_R1 | ![alt text](fastqc_plots/6_R1_per_base_qual.png) | ![alt text](qual_plots/6_R1_dist.png) |
| 6_R2 | ![alt text](fastqc_plots/6_R2_per_base_qual.png) | ![alt text](qual_plots/6_R2_dist.png) |

Here is a summary of runtime stats:
| Method | Runtime | Number of Cores | CPU Usage |
|---|---|---|---|
| FastQC | 2 minutes 49 seconds | 4 | 99% |
| qual_dist.py | 20 minutes 8 seconds | 8 | 99% |

The quality score distribution plots look very similar between the 2 methods. One obvious difference between the 2 plots is the the FastQC generated plot it gives more information about the distribution of quality scores for each nucleotide position than my plot. The red line shows the median, the yellow box shows the IQR, the upper and lower wisker give the 10% and 90% points, and the blue line gives the mean [1]. My plot only shows the mean for a given position, so I can compare the blue "mean" line with my plots. The lines for means appear to follow the same pattern, spiking and dipping in the same nucleotide positions, the dips and spikes are just less dramatic in the FastQC plots because the y-axis ranges from 0 to 30 whereas mine ranges from 30 to 40. The FastQC plots also show some guidelines with the red, orange, and green color zones to show what should be considered as a cutoff for 'good' quality, whereas my graph does not give any information about a cutoff. The FastQC plots are binned into 2 nucleotide groups, whereas my plots show a continuous distribution. This binning could potentially cause a loss of data. Finally, for my 4th file specifically, ```6_2D_mbnl_S5_L008_R2_001.fastq.gz```, the FastQC plot brings attention to some lower quality data: a dip in lower quartile/10% quality scores around bp 30-40, as well as a dip in lower quartile/10% quality at the end of the read from base pair 55 onwards. This means that there are some reads in the bottom 25% that are at or below the green cutoff of 28 in large sections of the read. My plot only shows the impact of these lower quality reads on the mean, which causes a small dip from about 38 to 36. I wouldn't think to throw out any reads based on my plot, but after seeing the FastQC  plot, I believe there are some lower quality reads in the ```6_2D_mbnl_S5_L008_R2_001.fastq.gz``` file that should be thrown out.

In terms of runtime, in addition to providing a more in depth analysis of the quality score distribution, FastQC also ran through all 4 files 10 times faster than my python script did (see above runtime comparison table). Additionally, the 2 methods used the same % of CPU but FastQC did so with half the number of cores as my script. Overall, FastQC provides a better analysis of quality distribution and is faster and more efficient.

Some reasons that FastQC might be so much better in performance and faster than my python script is FastQC is a java application [2], and java code runs faster than python code because java does not need go through an interpretation step like python does. Additionally, FastQC was first released in 2010 so has had many years and likely teams of people to optimize it, whereas my code was developed a couple weeks ago by 1 grad student who is new to bioinfomatics. Therefore, it makes sense why FastQC is better in every way.

source: [1] https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html
[2] https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt#:~:text=Installing%20FastQC%20%2D%2D%2D%2D%2D,Runtime%20Environment%20(JRE)%20installed.





Comment on the overall data quality of your two libraries. Go beyond per-base qscore distributions. Make and justify a recommendation on whether these data are of high enough quality to use for further analysis.








Questions:
* is my fastqc commend missing any parameters?
* did you do a quality cutoff to make these demuxed files?


### Part 2 - Adaptor trimming comparison
Goal: Perform adaptor trimming with existing tools cutadapt and Trimmomatic. 

First, I will install cutadapt and Trimmomatic into my QAA conda environment:
```
$ mamba install cutadapt=4.9
$ cutadapt --version
4.9
$ mamba install trimmomatic=0.39
$ trimmomatic -version
0.39
```
Now I will use cutadapt to trim adaptor sequences from my files.

Based on some research it appears that Illumina TruSeq uses the following adaptors:

Read 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
Read 2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

source: https://support-docs.illumina.com/SHARE/AdapterSequences/Content/SHARE/AdapterSeq/TruSeq/UDIndexes.htm

Sanity check: I will use bash commands to find the adaptor sequences.
```
zcat /projects/bgmp/shared/2017_sequencing/demultiplexed/16_3D_mbnl_S12_L008_R1_001.fastq.gz | grep -c "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" | less -p
zcat /projects/bgmp/shared/2017_sequencing/demultiplexed/16_3D_mbnl_S12_L008_R2_001.fastq.gz | grep "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
```
Sanity check: Use your Unix skills to search for the adapter sequences in your datasets and confirm the expected sequence orientations. Report the commands you used, the reasoning behind them, and how you confirmed the adapter sequences.

only get adaptors in seq if insert is shorter than read len. If insert < 101, you will get adaptor seq on 3' end .






I will use the following cutadapt command template:
```
cutadapt -a ADAPTER_FWD -A ADAPTER_REV -o out.1.fastq -p out.2.fastq reads.1.fastq reads.2.fastq
```
-a is the R1 adaptor sequence. -A is the R2 adaptor sequence. -o is the R1 output file. -p is the R2 output file. To run this on all my files, I put a command for each file pair into an sbatch script called ```run_cutadapt.sh``` and ran it.
```
$ sbatch run_cutadapt.sh 
Submitted batch job 15883794
```
The err and stout for this run is located in my ```Logs``` directory. The output files are located in the ```cutadapt_out``` directory.

Pair 16 Run Summary:
* 50 seconds to run, 6.166 µs/read
* used 98% CPU and 1 core
* Total read pairs processed: 8,235,197
* Read 1 with adapter: 1,002,983 (12.2%)
* Read 2 with adapter: 1,069,893 (13.0%)
* Pairs written (passing filters): 8,235,197 (100.0%)
* Total base pairs written (filtered):  1,627,540,628 bp (97.8%)
    * Read 1:   813,886,317/831,754,897 bp
    * Read 2:   813,654,311/831,754,897 bp
* Apdator 1 trimmed 1002983 times
* Adaptor 2 trimmed 1069893 times

Pair 6 Run Summary:
* 54 seconds to run, 4.904 µs/read
* used 98% CPU and 1 core
* Total read pairs processed: 11,028,244
* Read 1 with adapter: 416,045 (3.8%)
* Read 2 with adapter: 502,045 (4.6%)
* Pairs written (passing filters): 11,028,244 (100.0%)
* Total base pairs written (filtered):  2,220,049,445 bp (99.7%)
    * Read 1:   1,110,205,328/1,113,852,644 bp
    * Read 2:   1,109,844,117/1,113,852,644 bp
* Apdator 1 trimmed 416045 times
* Adaptor 2 trimmed 502045 times
	
---
### 8/26/24
### Part 2 cont
Now I will use Trimmomatic to quality trim my reads. I am going to write an sbatch script to run trimmomatic called ```run_trimmomatic.sh```.

I will use the following trimmomatic command template:
```
trimmomatic PE \
    <R1_file> <R2_file> \
    R1_paired_filtered.fastq.gz \
    R1_unpaired_filtered.fastq.gz \
    R2_paired_filtered.fastq.gz \
    R2_unpaired_filtered.fastq.gz \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:5:15 \
    MINLEN:35
```
* PE specifies parired end files
* LEADING is the threshold to remove low quality bases from the beginning of a read
* TRAILING is the threshold to remove low quality bases from the end of a read
* SLIDINGWINDOW specific sliding window trimming, where is the average quality within a window of 5 nucleotides is below 15 it cuts them out
* MINLEN is the length cutoff to remove reads that are short after cutting from the file

I ran this script with the following command:
```
$ sbatch run_trimmomatic.sh 
Submitted batch job 15884020
```

The log output is in the ```Logs``` directory. The outputted files are in the ```trim_out``` directory.

16 Summary:
* used 4 threads and 211% CPU
* took 3 minutes 10 seconds
* Input Read Pairs: 8235197
* Both Surviving: 8014158 (97.32%)
* Forward Only Surviving: 189384 (2.30%)
* Reverse Only Surviving: 6913 (0.08%)
* Dropped: 24742 (0.30%)

6 Summary
* used 4 threads and 200% CPU
* took 4 minutes 23 seconds
* Input Read Pairs: 11028244
* Both Surviving: 10461303 (94.86%)
* Forward Only Surviving: 552941 (5.01%)
* Reverse Only Surviving:  7434 (0.07%)
* Dropped: 6566 (0.06%)

It makes sense that there was a higher percent of only foward surviving reads in the 6 files because the R2 file was pretty low quality in my intial quality analysis.

Now, I will write a python script to plot the read length distributions for R1 and R2 reads for each file. It is called ```read_len_dist.py```.

To run the script, I needed to install ```matplotlib``` into my ```QAA``` environment
```
$ conda activate QAA
$ conda install matplotlib
```

I made an sbatch script called ```run_readlen_dist.sh``` to run my script on all the trimmomatic output files.
```
$ sbatch run_readlen_dist.sh 
Submitted batch job 15884140
```





 You can produce 2 different plots for your 2 different RNA-seq samples. There are a number of ways you could possibly do this. One useful thing your plot should show, for example, is whether R1s are trimmed more extensively than R2s, or vice versa. Comment on whether you expect R1s and R2s to be adapter-trimmed at different rates and why.

Questions:
* should we include the unpaired reads in our length distributions?






CHALLENGE - Run FastQC on your trimmed data. Comment on differences you observe between the trimmed and untrimmed data. Include any figures needed to support your conclusions.


,  and to demonstrate your ability to summarize other important information about this RNA-Seq data set in a high-level report. That is, you should create a cohesive, well written report for your "PI" about what you've learned about/from your data.