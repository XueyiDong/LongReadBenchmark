#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --output=ont_sub_cupcake_collapse.out
#SBATCH --mail-user=du.m@wehi.edu.au
#SBATCH --mail-type=END,FAIL

# Author: Mei Du
# Using cDNA_Cupcake version 25.2.0 (https://github.com/Magdoll/cDNA_Cupcake)
# Run on the downsampled dataset

SQANTI=/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/SQANTI3
CUPCAKE=/stornext/General/data/user_managed/grpu_mritchie_0/Mei/long_read_benchmark/SQANTI3/cDNA_Cupcake_dev

source activate $SQANTI/SQANTI3_env
export PYTHONPATH=$PYTHONPATH:$CUPCAKE/sequence
export PYTHONPATH=$PYTHONPATH:$CUPCAKE/rarefaction

collapse_isoforms_by_sam.py \
-b /stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/ONT_aligned/subsample_aligned/merged_ont_sub_10_sorted.bam \
--cpus 8 \
-o merged_ont_sub

