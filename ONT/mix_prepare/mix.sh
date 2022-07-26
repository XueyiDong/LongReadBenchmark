#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

cd /stornext/General/data/user_managed/grpu_mritchie_0/XueyiDong/long_read_benchmark/ONT/mix_fq

cat downsample/$1.fq.gz downsample/$2.fq.gz > $3.fq.gz