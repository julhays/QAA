#!/usr/bin/env python

# Author: Jules Hays

# PS8 Part 2: Parsing SAM file

import argparse

#define arg inputs, defaults are for the lane1_NoIndex_L001_R1_003.fastq file
def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Establish criteria for number of mapped and unmapped reads from star align output")
    parser.add_argument("-f", "--file_path", help="Specify the .sam star align output file", type=str, required=True)
    parser.add_argument("-o", "--out_num", help="Specify output filename prefix", type=str, required=True)
    return parser.parse_args()

#call get_args to create args object
args = get_args()

#set global variables and assign them to the user inputted values at the function call
file: str = args.file_path
num: str = args.out_num

#hard coded to test
#file = 'Danio_Aligned.out.sam'

#initialize counter variables to store counts
mapped_reads = 0
unmapped_reads = 0

#parse through the file
with open (file, 'rt') as fh:
    for line in fh:
        line = line.strip('\n')
        #skip passed header lines
        if line[0] == '@':
            continue
        #split each entry of the line into a list format
        line_sep = line.split('\t')

        #grab the bitwise flag in column 2
        flag = int(line_sep[1])

        #check to see if its a secondary alignment, if true don't do anything with it
        if((flag & 256) != 256):
            #check to see if read is unmapped - true (not 4) means mapped and false means unmapped
            if((flag & 4) != 4):
                mapped_reads += 1
            else:
                unmapped_reads += 1



#print the results
print(f'Results for file {num}')
print("Number of mapped reads:", mapped_reads)
print("Number of unmapped reads:", unmapped_reads)

with open(f'mapped_counts/{num}_mapped_reads.txt', 'wt') as out:
    out.write(f'Number of mapped reads: {mapped_reads}\n')
    out.write(f'Number of unmapped reads: {unmapped_reads}')

exit()
            

            
