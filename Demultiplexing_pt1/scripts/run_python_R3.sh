#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --job-name=demu_R3
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




file3="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"






/usr/bin/time -v python ./c_count.py -f $file3 -o 1294_S1_L008_R3_001_hist
