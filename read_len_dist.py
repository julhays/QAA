#!/usr/bin/env python

# Author: Jules Hays

# QAA Part 2: Read Length Distribution after Quality Trimming

import argparse
import matplotlib.pyplot as plt
import gzip

#define arg inputs, defaults are for the lane1_NoIndex_L001_R1_003.fastq file
def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Establish criteria for generating read length distribution plot")
    parser.add_argument("-p1", "--R1_paired", help="Specify the path of the R1 paired trimmomatic output", type=str, required = True)
    parser.add_argument("-u1", "--R1_unpaired", help="Specify the path of the R1 paired trimmomatic output", type=str, required = True)
    parser.add_argument("-p2", "--R2_paired", help="Specify the path of the R1 paired trimmomatic output", type=str, required=True)
    parser.add_argument("-u2", "--R2_unpaired", help="Specify the path of the R1 paired trimmomatic output", type=str, required=True)
    parser.add_argument("-o", "--output_name", help="Specify the name of the output graph with no .png extension", type=str, required=True)
    return parser.parse_args()

#call get_args to create args object
args = get_args()

#set global variables and assign them to the user inputted values at the function call
R1_paired: str = args.R1_paired
R1_unpaired: str = args.R1_unpaired
R2_paired: str = args.R2_paired
R2_unpaired: str = args.R2_unpaired
output: str = args.output_name

#hard code file paths for testing
# R1_paired = 'trim_out/16_R1_paired_filtered.fastq.gz'
# R1_unpaired = 'trim_out/16_R1_unpaired_filtered.fastq.gz'
# R2_paired = 'trim_out/16_R2_paired_filtered.fastq.gz'
# R2_unpaired = 'trim_out/16_R2_unpaired_filtered.fastq.gz'
# output = '16'

#functions
def read_len_dist(file: str, len_list: list):  #no return
    '''Takes in a fastq file and a list to store lengths, 
    appends lengths of each line onto the list. Returns nothing'''
    with gzip.open(file, "rt") as fh:
        counter = 0 #counter to track line number
        #loop through input fastq file
        for line in fh:
            counter+=1
            line = line.strip('\n')
            if counter%2 == 0:              #pulls only the sequence line
                #append the length of the sequence to the respective list
                len_list.append(len(line))
    
    return

#initialize empty lists to store read lengths in the file
r1_lens = []
r2_lens = []

#add the distribution of read lengths for all files to respective lists
r1p_lens = read_len_dist(R1_paired, r1_lens)
r1u_lens = read_len_dist(R1_unpaired, r1_lens)
r2p_lens = read_len_dist(R2_paired, r2_lens)
r2u_lens = read_len_dist(R2_unpaired, r2_lens)

#plot the distribution of read lengths for R1 and R2
input = [r1_lens, r2_lens]
plt.hist(input, histtype='bar')
plt.xlabel("Read length (nts)")
plt.ylabel("Frequency")
plt.title(f'Read Length Distribution for Trimmed File Pair {output}')
plt.legend(('R1', 'R2'))
plt.savefig(f'read_len_plots/read_len_dist_{output}.png')


