#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL


outdir=/vast/scratch/users/dong.x/long_read_benchmark/pilot_ONT/bam_tx
unixhome=/wehisan/home/allstaff/d/dong.x
ref=/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/HumanSequins/gencode.v33.sequins.transcripts.fa
export PATH=$PATH:/wehisan/home/allstaff/d/dong.x/Programs/minimap2

module load samtools/1.7
# 1: fq dir 2: fq file name 
minimap2 -ax map-ont -t 8 --sam-hit-only $ref $1/$2.fq.gz | samtools view -b | samtools sort > $outdir/$2.sorted.bam
