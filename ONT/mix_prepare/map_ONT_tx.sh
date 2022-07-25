#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=300G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

cd /stornext/General/data/user_managed/grpu_mritchie_0/XueyiDong/long_read_benchmark/ONT/mix_fq
mkdir -p bam_tx
unixhome=/wehisan/home/allstaff/d/dong.x
# fq=/stornext/Projects/promethion/promethion_access/lab_ritchie/transcr_bench/long_term/transcr_bench/guppy_4_rebasecall_Oct20

ref=/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/HumanSequins/gencode.v33.sequins.transcripts.fa
export PATH=$PATH:/wehisan/home/allstaff/d/dong.x/Programs/minimap2

module load samtools/1.7

# for sample in barcode{05..06}
# do minimap2 -ax map-ont -t 8 --sam-hit-only $ref $fq/$sample.fq.gz | samtools view -b | samtools sort > $sample.bam
# samtools index $sample.sorted.bam $sample.bai
# done

minimap2 -ax map-ont -t 8 --sam-hit-only $ref ./$1.fq.gz | samtools view -b | samtools sort > bam_tx/$1.bam
