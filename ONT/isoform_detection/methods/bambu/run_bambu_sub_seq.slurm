#!/bin/bash
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=30G
#SBATCH --job-name=bambu
#SBATCH --output=ONT_sub_bambu.out
#SBATCH --mail-user=du.m@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

# Author: Mei Du
# Using bambu version 1.0.3 (https://github.com/GoekeLab/bambu)
# Run on the downsampled dataset with the modified GTF

module load R/4.1.1
Rscript /stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/bambu/scripts/run_bambu_sub_seq.R