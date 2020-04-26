#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --job-name=demu_R4
#SBATCH --output=slurm-%j-%x.out

#SBATCH --time=0-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

conda deactivate
conda deactivate
conda deactivate
conda deactivate
conda activate bgmp_py3




file4="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"






/usr/bin/time -v python ./c_count.py -f $file4 -o 1294_S1_L008_R4_001_hist
