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


for f in *.fastq.gz ; do
  mv -- "$f" "${f%.fastq.gz}.fastq"
  done
/usr/bin/time -v gzip *.fastq
