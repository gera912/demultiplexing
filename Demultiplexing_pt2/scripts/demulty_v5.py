#!/usr/bin/env python

# Program Header
# Course: Bi624
# Name:   Gerardo Perez
# Description: Our goal is to look through a lane of sequencing library
# preps and sort pooled experiments of a lane of Illumina sequencing data.
# The script will ask to input 4 fastq files to be read. The output will be a total of
# 52 fastq.gz files: R1 matching pair index file, R2 matching pair index file, R1
# undertermine list file, R2 undertermine list file, R1 index hopped list file and R2 index hopped list file.
# In addition, there will be an output that will have percent frequencies.
#
#
# Demultiplexing part 2
#
#
# Development Environment: Atom 1.38.2
# Version: Python 3.7.3
# Date:  11/05/2019
#################################################

# Imports module
import matplotlib.pyplot as plt
import argparse
import gzip
import numpy as np
import statistics


# Creates an arguement passer
parser = argparse.ArgumentParser(description="A program for base distribution")


# Adds arguemets by calling the arguement passer
parser.add_argument("-f1", "--filename1", help="specify the filename1", required=True)
#args = parser.parse_args()

# Adds arguemets by calling the arguement passer
parser.add_argument("-f2", "--filename2", help="specify the filename2", required=True)
#args = parser.parse_args()

# Adds arguemets by calling the arguement passer
parser.add_argument("-f3", "--filename3", help="specify the filename3", required=True)
#args = parser.parse_args()

# Adds arguemets by calling the arguement passer
parser.add_argument("-f4", "--filename4", help="specify the filename4", required=True)
#args = parser.parse_args()

# Parses arguemets through using the parse args method.
args = parser.parse_args()

# creates variables for args parse
file1=args.filename1

file2=args.filename2

file3=args.filename3

file4=args.filename4

# creates variables to store files that are open to be written
U_R1=open("Undetermined_list_R1.fastq.gz"  , "w")
U_R2=open("Undetermined_list_R2.fastq.gz"  , "w")
I_R1=open("Index_Hopped.fastq_R1.fastq.gz"  , "w")
I_R2=open("Index_Hopped.fastq_R2.fastq.gz"  , "w")

# creates a variable to store the index file name for barcodes
barcodes_list="indexes.txt"

# Craetes a function to get the reverse complement of a sequence
def rev_comp(seq):
    """Takes a sequence and returns the reverse the complement"""
    valid_bases = ['A', 'T', 'G', 'C']

    #Creates a dictionary for nucleotide complementary bases.
    comp_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    bases=[]

    # A for loop that checks the characters in the sequece string input.
    for base in seq:

        # If statement that checks if the character is in the list of valid bases.
        if base in valid_bases:

            # If true then get the complementary base from dictionary of complementary bases.
            bases.append(comp_dict[base])

        # If false, then add base to the list of bases.
        else: bases.append(base)

    # returns the list as a string
    rev_seq=bases[::-1]

    # returns the list as a string
    return "".join(rev_seq)

# Creates a method that checks if the average socre of the bases in the biological sequence reads is less than the cutoff 30
# The reason for the 30 cutoff is due from looking at the average means scores histogram in the part 1 biological reads.
# This score cutoff will include almost all the average base pairs for each position and not throw out too much data.
def ave_base_qc(qc_line):
    """Takes a sequence and returns false if the average base in sequnece is less than 30"""

    # Creates an empty list to store values.
    values_list=[]

    # Creates an empty list to store means.
    mean_list=[]

    i = 0

    # Creates a for loop to check the character of the sequence string input
    for char in qc_line:

        #Calls a function to convert the character into a phred score and stores the result to a variable
        value=convert_phred(char)

        # Adds the converted phred value into a list
        values_list.append(value)

    # Calculates the average of the values in list
    mean_value=sum(values_list)/len(values_list)

    # If statenent that checks if mean value is less than 30
    if mean_value<30:

        # If true, then return false.
        return False




# Opens a file to read barcode list than stores this to a variable
f=open(barcodes_list,"r")

# Stores the lines from barcode list to a variable.
lines=f.readlines()

# Creates an empty list for barcodes.
barcodes=[]
count=0

# A for loop that goes through each line.
for i in lines:

    # Strips the new line charcter, splits through each tab and gets the 4th column to store in a value
    x=i.strip().split('\t')[4]

    # Adds the value that was in the 4th column of each line into a list
    barcodes.append(x)

barcode_dict={}

# Removes the first index of the list
barcodes=barcodes[1:]


barcode_count={}

# A for loop that goes through each barcode of barcode list
for barcode in barcodes:

    # Creates the string for the barcode dual match to a variable.
    x=barcode+"-"+barcode

    # Sets the barcode dual match to dictionary with 0 as values
    barcode_count.setdefault(x,0)

    # Stores the string for the barcode dual match variable as a key
    # and opens R1 and R2, with barcode dual match names, to be written as values
    barcode_dict[x]=[open(x+"_R1.fastq.gz", "w"), open(x+"_R2.fastq.gz", "w")]




f.close()

# Creates variables
LN = 0
undetermine=0
index_hopped=0
index_paired=0
low_quality_barcode=0
low_quality_reads=0

list_files=[file1, file2, file3, file4]


# Opens a gz files through the use of argparse
r1=gzip.open(file1, "rt")
r2=gzip.open(file2, "rt")
r3=gzip.open(file3, "rt")
r4=gzip.open(file4, "rt")


# Creates a method that checks if a base score in barcode sequence reads is less than the cutoff 20
# The reason for the 20 cutoff is due to not throw out too much data.
def barcode_qc(qc_line):
    """Takes a sequence and returns false if quality control barcode base cutoff is less than 20"""

    for char in qc_line:

        value=convert_phred(char)

        if value < 20:

            return False


def convert_phred(letter):
    """Converts a single character into a phred score"""
    x = ord(letter) - 33
    return x

# A while loop use to read files
while True:

    # Creates lists to save lines from fasta file and remove them after every iteration
    file1_lines=[]
    file2_lines=[]
    file3_lines=[]
    file4_lines=[]

    # For loop that will iterate 4 times to store the 4 lines from file into a list.
    # breaks loop if empty.
    for i in range(4):
        f1_R1=r1.readline().strip()
        if f1_R1=="":
            break
        file1_lines.append(f1_R1)
    if f1_R1=="":
        break

    for i in range(4):
        f1_R2=r2.readline().strip()
        if f1_R2=="":
            break
        file2_lines.append(f1_R2)
    if f1_R2=="":
        break

    for i in range(4):
        f1_R3=r3.readline().strip()
        if f1_R3=="":
            break
        file3_lines.append(f1_R3)
    if f1_R3=="":
        break


    for i in range(4):
        f1_R4=r4.readline().strip()
        if f1_R4=="":
            break
        file4_lines.append(f1_R4)
    if f1_R4=="":
        break

    # For loop that will iterate 4 times.
    for i in range(1,len(file1_lines), 4):

            # If statement that checks if there is Ns in barcode sequence in R2 and R3 files,
            # if true then writes to R1 and R2 undetermine files
            if  "N" in file2_lines[i] or "N" in file3_lines[i]:

                U_R1.write("{0}\t{1}-{2}\n{3}\n{4}\n{5}\n".format(file1_lines[i-1],file2_lines[i], file3_lines[i],file1_lines[i],file1_lines[i+1],file1_lines[i+2]))
                U_R2.write("{0}\t{1}-{2}\n{3}\n{4}\n{5}\n".format(file4_lines[i-1],file2_lines[i], file3_lines[i],file4_lines[i],file4_lines[i+1],file4_lines[i+2]))

                undetermine=undetermine+1

            # else if statement that checks the barcode quality control barcode cutoff 20 in R2 and R3 files
            # if one is less then writes to R1 and R2 undetermine files.
            elif (barcode_qc(file2_lines[i+2])== False) or (barcode_qc(file3_lines[i+2])== False):
                U_R1.write("{0}\t{1}-{2}\n{3}\n{4}\n{5}\n".format(file1_lines[i-1],file2_lines[i], file3_lines[i],file1_lines[i],file1_lines[i+1],file1_lines[i+2]))
                U_R2.write("{0}\t{1}-{2}\n{3}\n{4}\n{5}\n".format(file4_lines[i-1],file2_lines[i], file3_lines[i],file4_lines[i],file4_lines[i+1],file4_lines[i+2]))

                low_quality_barcode=low_quality_barcode+1

            # else if statement that checks the average base quality control barcode cutoff 30 in R1 and R4 files
            # if one is less then writes to R1 and R2 undetermine files.
            elif (ave_base_qc(file1_lines[i+2])==False) or (ave_base_qc(file4_lines[i+2])==False):
                U_R1.write("{0}\t{1}-{2}\n{3}\n{4}\n{5}\n".format(file1_lines[i-1],file2_lines[i], file3_lines[i],file1_lines[i],file1_lines[i+1],file1_lines[i+2]))
                U_R2.write("{0}\t{1}-{2}\n{3}\n{4}\n{5}\n".format(file4_lines[i-1],file2_lines[i], file3_lines[i],file4_lines[i],file4_lines[i+1],file4_lines[i+2]))

                low_quality_reads=low_quality_reads+1

            # else if statement that checks if barcodes in R2 or reverse complement in R3 files are not in barcodes list.
            # if true then writes to R1 and R2 undetermine files.
            elif (file2_lines[i] not in barcodes) or (rev_comp(file3_lines[i]) not in barcodes):
                U_R1.write("{0}\t{1}-{2}\n{3}\n{4}\n{5}\n".format(file1_lines[i-1],file2_lines[i], file3_lines[i],file1_lines[i],file1_lines[i+1],file1_lines[i+2]))
                U_R2.write("{0}\t{1}-{2}\n{3}\n{4}\n{5}\n".format(file4_lines[i-1],file2_lines[i], file3_lines[i],file4_lines[i],file4_lines[i+1],file4_lines[i+2]))

                undetermine=undetermine+1

            # else if statement that checks if barcodes in R2 equals the reverse complement barcodes in R3 files and R2 barcodes are in barcode list.
            # if true then writes to R1 and R2 index paired files.
            elif (file2_lines[i]==rev_comp(file3_lines[i])) and (file2_lines[i] in barcodes):

                match=file2_lines[i]+"-"+rev_comp(file3_lines[i])
                barcode_dict[match]=open(match+"_R1.fastq.gz", "a").write("{0}\t{1}-{2}\n{3}\n{4}\n{5}\n".format(file1_lines[i-1],file2_lines[i], file3_lines[i],file1_lines[i],file1_lines[i+1],file1_lines[i+2]))
                barcode_dict[match]=open(match+"_R2.fastq.gz", "a").write("{0}\t{1}-{2}\n{3}\n{4}\n{5}\n".format(file1_lines[i-1],file2_lines[i], file3_lines[i],file1_lines[i],file1_lines[i+1],file1_lines[i+2]))


                index_paired=index_paired+1
                barcode_count[match]=barcode_count[match]+1

            # else statement that will write to R1 and R2 index hopped files.
            else:

                I_R1.write("{0}\t{1}-{2}\n{3}\n{4}\n{5}\n".format(file1_lines[i-1],file2_lines[i], file3_lines[i],file1_lines[i],file1_lines[i+1],file1_lines[i+2]))
                I_R2.write("{0}\t{1}-{2}\n{3}\n{4}\n{5}\n".format(file4_lines[i-1],file2_lines[i], file3_lines[i],file4_lines[i],file4_lines[i+1],file4_lines[i+2]))

                index_hopped=index_hopped+1
            count=count+1

# Closes files.
open("Undetermined_list_R1.fastq.gz"  , "w").close()
open("Undetermined_list_R2.fastq.gz"  , "w").close()
open("Index_Hopped.fastq_R1.fastq.gz"  , "w").close()
open("Index_Hopped.fastq_R2.fastq.gz"  , "w").close()

# For loop that will print the dual match barcodes with the percentage frequency.
for key, value in barcode_count.items():
    print(key, round((value/count)*100,1),"%")

# Print statements that will print the percentage frequencies
print(sum(barcode_count.values()))
print("Reads"+":",count)
print("undetermine"+":",round((undetermine/count)*100,1), "%")
print("index_hopped"+":",round((index_hopped/count)*100,1), "%")
print("index_paired"+":",round((index_paired/count)*100,1), "%")
print("low_quality_barcode"+":",round((low_quality_barcode/count)*100,1), "%")
print("low_quality_reads"+":",round((low_quality_reads/count)*100,1),"%")


# Program Output (Commented out)
"""
GTAGCGTA-GTAGCGTA 1.8 %
CGATCGAT-CGATCGAT 1.3 %
GATCAAGG-GATCAAGG 1.4 %
AACAGCGA-AACAGCGA 2.0 %
TAGCCATG-TAGCCATG 2.3 %
CGGTAATC-CGGTAATC 0.9 %
CTCTGGAT-CTCTGGAT 7.7 %
TACCGGAT-TACCGGAT 16.0 %
CTAGCTCA-CTAGCTCA 4.1 %
CACTTCAC-CACTTCAC 0.9 %
GCTACTCT-GCTACTCT 1.5 %
ACGATCAG-ACGATCAG 1.9 %
TATGGCAC-TATGGCAC 2.5 %
TGTTCCGT-TGTTCCGT 3.6 %
GTCCTAAG-GTCCTAAG 2.0 %
TCGACAAG-TCGACAAG 0.9 %
TCTTCGAC-TCTTCGAC 9.3 %
ATCATGCG-ATCATGCG 2.3 %
ATCGTGGT-ATCGTGGT 1.5 %
TCGAGAGT-TCGAGAGT 2.5 %
TCGGATTC-TCGGATTC 1.0 %
GATCTTGC-GATCTTGC 0.8 %
AGAGTCCA-AGAGTCCA 2.5 %
AGGATAGC-AGGATAGC 1.9 %
263427112
Reads: 363246735
undetermine: 2.9 %
index_hopped: 0.1 %
index_paired: 72.5 %
low_quality_barcode: 21.6 %
low_quality_reads: 2.9 %
	Command being timed: "python demulty_v5.py -f1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -f2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -f3 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -f4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"
	User time (seconds): 23995.37
	System time (seconds): 6240.21
	Percent of CPU this job got: 64%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 12:57:27
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 161936
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 328180
	Voluntary context switches: 4680130
	Involuntary context switches: 874954
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0

"""
