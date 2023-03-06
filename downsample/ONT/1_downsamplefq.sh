#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --job-name=downsample
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=END,FAIL

# 1 sample 2 number of reads 3 output file name 4 input fastq directory



# unixhome=/wehisan/home/allstaff/d/dong.x
# FQ=/stornext/Projects/promethion/promethion_access/lab_ritchie/transcr_bench/long_term/transcr_bench/guppy_4_rebasecall_Oct20/barcode01-06

export PATH=$PATH:/wehisan/home/allstaff/d/dong.x/Programs/seqtk

seqtk sample -s100 $4/$1.fq.gz $2 | gzip > $3.fq.gz
