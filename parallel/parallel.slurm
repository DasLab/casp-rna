#!/bin/bash
#SBATCH -p owners,biochem
#SBATCH -t 24:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mail-type=FAIL,END
#################
#a file for job output, you can check job progress
#SBATCH --output=casp-rna_%j.out
#################
# a file for errors from the job
#SBATCH --error=casp-rna_%j.err
#################

target=$1
metric=$2

ml py-pandas/1.3.1_py39;  ml py-numpy/1.20.3_py39; ml viz; ml py-matplotlib/3.4.2_py39; ml py-scipy/1.6.3_py39
python run_parallel.py $target $metric