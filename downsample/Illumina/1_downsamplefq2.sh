#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --job-name=Illumina_downsample
#SBATCH --mail-type=END,FAIL

# 1 sample 2 number of reads 3 output file name 4 input fastq directory

# cd /stornext/General/data/user_managed/grpu_mritchie_0/XueyiDong/long_read_benchmark/illumina
# mkdir -p mix_fq/downsample

# unixhome=/wehisan/home/allstaff/d/dong.x
# FQ=/stornext/General/data/user_managed/grpu_mritchie_1/long_RNA_benchmark/Illumina_April20/data/raw/

export PATH=$PATH:/wehisan/home/allstaff/d/dong.x/Programs/seqtk

seqtk sample -s100 $4/$1.fq.gz $2 | gzip > $3.fq.gz
