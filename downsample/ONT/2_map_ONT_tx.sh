#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --job-name=map_ONT
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=END,FAIL

outdir=/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/downsample/ONT/bam
fqdir=/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/downsample/ONT/fastq
ref=/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/HumanSequins/gencode.v33.sequins.transcripts.fa
export PATH=$PATH:/wehisan/home/allstaff/d/dong.x/Programs/minimap2

module load samtools/1.7
module load stornext
# 1 fq file name 
minimap2 -ax map-ont -t 16 --sam-hit-only $ref $fqdir/$1.fq.gz | samtools view -b | samtools sort > $outdir/$1.bam
# truncate fastq after mapping
# snrmdiskcopy $fqdir/$1.fq.gz
