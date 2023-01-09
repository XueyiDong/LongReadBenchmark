#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

# cd /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/ONT

export PATH=/home/users/allstaff/dong.x/Programs/salmon/salmon-latest_linux_x86_64/bin:$PATH

module load samtools
TX=/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/HumanSequins/gencode.v33.sequins.transcripts.fa

# 1: out directory 2: input BAM file (without .bam)

# for SAMPLE in barcode{01..06}
# do 
# mkdir -p ./salmon_bs/$SAMPLE
samtools view -F 256 -b -@ 8 $2.bam > $2.primary.bam
salmon quant -t $TX -l A -a $2.primary.bam -o $1 --numBootstraps 100 -p 8
# done