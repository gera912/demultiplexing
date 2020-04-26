#!/usr/bin/env python

# Program Header
# Course: Bi622
# Name:   Gerardo Perez
# Description: Our goal with this assignment is to assess the overall quality
# of a lane of Illumina sequencing data. Our goal is to calculate the average
# quality score along each of the basepair of data we have. Then output is a
# histogram of average base pair means.
#
# Demultiplexing: c_count.py
#
#
# Development Environment: Atom 1.38.2
# Version: Python 3.7.3
# Date:  06/19/2019
#################################################

# Imports module
import matplotlib.pyplot as plt
import argparse
import gzip

#file="test_R1.txt.gz"
#file="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"



# Creates an arguement passer
parser = argparse.ArgumentParser(description="A program for base distribution")


# Adds arguemets by calling the arguement passer
parser.add_argument("-f", "--filename", help="specify the filename", required=True)
#args = parser.parse_args()

#Adds arguemets by calling the arguement passer
parser.add_argument("-o", "--outputname", help="specify the filename", required=True)

# Parses arguemets through using the parse args method.
args = parser.parse_args()

# creates variables for args parse
file=args.filename
output=args.outputname


# Creates a list
mean_scores = []

# Initiates a variable.
LN = 0

# opens and unzips file to read and stores the file as a variable.
with gzip.open(file,"rt") as fh:

    # A while loop.
    while True:

        # Reads the line and strips the new line character.
        # Stores the result to a variable. Moves to next line.
        fh.readline().strip()

        # updates the variable by an addition of one after each line has been read.
        LN = 1 + LN

        # Reads the line and strips the new line character.
        # Stores the result to a variable. Moves to next line.
        L2=fh.readline().strip()

        # Stores the lenghth of line 2 into a variable.
        total=len(L2)

        # Break loop.
        break

# Creates a method that takes an instance variable.
def convert_phred(letter):
    """Converts a single character into a phred score"""

    # Creates a variable that stores a difference
    # of -33 of the unicode point from the instance variable.
    x = ord(letter) - 33

    # returns a variable with the result.
    return x

# Creates an empty list
mean_scores = []

# A for loop to insert 101 elements.
for i in range(total):

    # Adds a float to the list.
    mean_scores.append(0.0)


    # Initiates a variable.
    LN = 0

# opens and unzips file to read and stores the file as a variable.
with gzip.open(file,"rt") as fh:

     # A for loop that goes through each line in the text file.
    for line in fh:

        # updates the variable by an addition of one after each loop.
        LN = LN + 1

        # If statemnt to get the fourth line from each record in text file.
        if LN%4 == 0:

            # A for loop that goes through each character in the line.
            for i in range(total):

                # Stores a character from line to a variable
                char = line[i]

                # calls the method with an instance ands stores the
                # result into a variable.
                value = convert_phred(char)

                #  stores and sums up values into a list
                mean_scores[i] = mean_scores[i] + value

# Stores the value of number of records to a varaible.
Num_rec = LN/4


# A for loop that goes through each index in the length of the list.
for i in range(len(mean_scores)):

    # Divides each index in the list by the number of records and stores
    # the result back to list.
    mean_scores[i] = (mean_scores[i])/(Num_rec)

# Prints the length of number of bases.
# print(len(mean_scores))

# Plots a bar chart
plt.bar(range(len(mean_scores)), mean_scores)

# Labels X axis
plt.xlabel('Base Pair')

# Labels Y axis
plt.ylabel('Mean Quality Score')

# outputs a png file
plt.savefig(output)
