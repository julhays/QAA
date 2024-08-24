# BI623 RNA-seq Quality Assessment Assignment - Lab Notebook
Jules Hays
## Due: 9/13/24

Python version: 3.12

Environments used: QAA (fastqc 0.12.1, cutadapt 4.9, trimmomatic 0.39)

bgmp_py312 (matplotlib 3.9.1), 
bgmp_velvet (Velvet 1.2.10)

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

The files are located at the followng paths:
```
/projects/bgmp/shared/2017_sequencing/demultiplexed/16_3D_mbnl_S12_L008_R1_001.fastq.gz
/projects/bgmp/shared/2017_sequencing/demultiplexed/16_3D_mbnl_S12_L008_R2_001.fastq.gz
/projects/bgmp/shared/2017_sequencing/demultiplexed/6_2D_mbnl_S5_L008_R1_001.fastq.gz
/projects/bgmp/shared/2017_sequencing/demultiplexed/6_2D_mbnl_S5_L008_R2_001.fastq.gz
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

| File name | size (M) | # lines | # records | read_len | Phred encoding |
|---|---|---|---|---|---|
| 16_3D_mbnl_S12_L008_R1_001.fastq.gz | 468 | 32,940,788 | 8,235,197 | 101 | Phred+33 |
| 16_3D_mbnl_S12_L008_R2_001.fastq.gz | 452 | 32,940,788 | 8,235,197 | 101 | Phred+33 |
| 6_2D_mbnl_S5_L008_R1_001.fastq.gz | 571 | 44,112,976  | 11,028,244 | 101 | Phred+33 |
| 6_2D_mbnl_S5_L008_R2_001.fastq.gz | 653 | 44,112,976  | 11,028,244 | 101 | Phred+33 |

The files look good to continue with quality analysis.


### Part 1 - Read quality score ditributions
Goal: Install and use FastQC for quality assessment of demultiplexed RNA-seq files

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
| 16_3D_mbnl_S12_L008_R1_001 | ![alt text](notebook_images/16_R1_per_base_qual.png) | ![alt text](notebook_images/16_R1_per_base_n.png) |
| 16_3D_mbnl_S12_L008_R2_001 | ![alt text](notebook_images/16_R2_per_base_qual.png) | ![alt text](notebook_images/16_R2_per_base_n.png) |
| 6_2D_mbnl_S5_L008_R1_001 | ![alt text](notebook_images/6_R1_per_base_qual.png) | ![alt text](notebook_images/6_R1_per_base_n.png) |
| 6_2D_mbnl_S5_L008_R2_001 | ![alt text](notebook_images/6_R2_per_base_qual.png) | ![alt text](notebook_images/6_R2_per_base_n.png) |


Questions:
* is my fastqc commend missing any parameters?



### Part 2 - Adaptor trimming comparison

GOal: Perform adaptor trimming with existing tools cutadapt and Trimmomatic. 
, compare the quality assessments to those from your own software, and to demonstrate your ability to summarize other important information about this RNA-Seq data set in a high-level report. That is, you should create a cohesive, well written report for your "PI" about what you've learned about/from your data.

First, I will install cutadapt and Trimmomatic into my QAA conda environment:
```
$ mamba install cutadapt=4.9
$ cutadapt --version
4.9
$ mamba install trimmomatic=0.39
$ trimmomatic -version
0.39
```
