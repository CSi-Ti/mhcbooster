#!/bin/bash

#SBATCH --job-name=mhcv
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --partition=day
#SBATCH --cpus-per-task=16
#module --force purge
#source ~/miniconda3/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"

conda activate mhcvalidator

python /home/rw762/workspace/mhc-validator/run_mhcvalidator.py
