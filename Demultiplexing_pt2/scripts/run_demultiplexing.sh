#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --job-name=Demulti_trial
#SBATCH --output=slurm-%j-%x.out

#SBATCH --time=0-24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=7

conda deactivate
conda deactivate
conda deactivate
conda deactivate
conda activate bgmp_py3

file1="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
file2="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
file3="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
file4="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"


/usr/bin/time -v python demulty_v5.py -f1 $file1 -f2 $file2  -f3 $file3 -f4 $file4
