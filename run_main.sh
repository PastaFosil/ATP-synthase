#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --account=def-dsivak
#SBATCH --ntasks=10               # number of MPI processes
#SBATCH --mem-per-cpu=1G      # memory; default unit is megabytes
module load python 
module load scipy-stack
python3 "/home/juanch0/projects/def-dsivak/juanch0/ATP_ratio/main.py" 8.0
