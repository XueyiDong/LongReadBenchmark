#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

# 1 sample 2 number of reads 3 output file name

cd /stornext/General/data/user_managed/grpu_mritchie_0/XueyiDong/long_read_benchmark/ONT
mkdir -p mix_fq/downsample

unixhome=/wehisan/home/allstaff/d/dong.x
FQ=/stornext/Projects/promethion/promethion_access/lab_ritchie/transcr_bench/long_term/transcr_bench/guppy_4_rebasecall_Oct20/barcode01-06

export PATH=$PATH:/wehisan/home/allstaff/d/dong.x/Programs/seqtk

seqtk sample -s100 $FQ/$1.fq.gz $2 | gzip > mix_fq/downsample/$3.fq.gz

# for SAMPLE in barcode{01..06}
# do seqtk sample -s100 $SAMPLE.fq.gz 0.2 | gzip > mix_fq/downsample/$SAMPLE.fq.gz
# done