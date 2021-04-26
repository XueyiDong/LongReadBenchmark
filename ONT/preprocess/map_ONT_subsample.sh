#!/bin/bash
#PBS -q submit
#PBS -l nodes=1:ppn=8,mem=80gb,walltime=480:00:00
#PBS -M dong.x@wehi.edu.au
#PBS -m abe
#PBS -j oe

cd /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/ONT/bam_subsample
unixhome=/wehisan/home/allstaff/d/dong.x
# fq=/stornext/Projects/promethion/promethion_access/lab_ritchie/transcr_bench/long_term/transcr_bench/fq_merged
fq=/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/subsample/ONT010
genome=$unixhome/annotation/HumanSequins/GrCh38_sequins.fa
bed=$unixhome/annotation/HumanSequins/gencode.v33.sequins.junction.bed
export PATH=$PATH:/wehisan/home/allstaff/d/dong.x/Programs/minimap2

module load samtools/1.7

for sample in barcode{01..06}
do minimap2 -ax splice -uf -k14 --junc-bed  $bed $genome $fq/$sample.fq.gz | samtools view -b | samtools sort > $sample.sorted.bam
samtools index $sample.sorted.bam $sample.sorted.bai
done