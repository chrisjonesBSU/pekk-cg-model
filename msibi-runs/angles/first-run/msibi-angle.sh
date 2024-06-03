#!/bin/bash
#SBATCH --gres gpu:1
#SBATCH --output=log.o
#SBATCH --error=log.e
#SBATCH --ntasks=1
source /home/chrisjones4/.bashrc
conda activate msibi 
python angles-msibi.py
