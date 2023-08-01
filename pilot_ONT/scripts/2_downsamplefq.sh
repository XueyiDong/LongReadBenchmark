#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --job-name=downsample
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

# 1 sample (barcodeXX) 2 number of reads 3 output file name

cd /vast/scratch/users/dong.x/long_read_benchmark/pilot_ONT/
mkdir -p fq_subsample


fq=/stornext/Projects/promethion/promethion_access/lab_ritchie/Xueyi_cellLine_mixture/long_term/fastq/H1_1-HC_5/20221107_1451_2D_PAK15503_0f50058f/guppy6.2.1_sup_prom

export PATH=$PATH:/wehisan/home/allstaff/d/dong.x/Programs/seqtk

seqtk sample -s100 $fq/$1\_sup_q10_pass.fq.gz $2 | gzip > fq_subsample/$3.fq.gz