#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

cd /stornext/General/data/user_managed/grpu_mritchie_0/XueyiDong/long_read_benchmark/illumina/mix_fq

cat downsample/$1.fq.gz downsample/$2.fq.gz > $3.fq.gz

# cat R1-000-sub25.fastq R1-100-sub75.fastq > R1-075sub.fastq
# cat R1-000-sub50.fastq R1-100-sub50.fastq > R1-050sub.fastq
# cat R1-000-sub75.fastq R1-100-sub25.fastq > R1-025sub.fastq
# cat R2-000-sub25.fastq R2-100-sub75.fastq > R2-075sub.fastq
# cat R2-000-sub50.fastq R2-100-sub50.fastq > R2-050sub.fastq
# cat R2-000-sub75.fastq R2-100-sub25.fastq > R2-025sub.fastq
# cat R3-000-sub25.fastq R3-100-sub75.fastq > R3-075sub.fastq
# cat R3-000-sub50.fastq R3-100-sub50.fastq > R3-050sub.fastq
# cat R3-000-sub75.fastq R3-100-sub25.fastq > R3-025sub.fastq