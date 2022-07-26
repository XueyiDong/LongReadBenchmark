#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

cd /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/illumina/
mkdir -p salmon_bs


export PATH=/home/users/allstaff/dong.x/Programs/salmon/salmon-latest_linux_x86_64/bin:$PATH

TX=/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/HumanSequins/gencode.v33.sequins.transcripts.fa
FQ=/stornext/General/data/user_managed/grpu_mritchie_0/XueyiDong/long_read_benchmark/illumina/mix_fq

mkdir -p ./salmon_bs/$1
salmon quant -i /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/HumanSequins/salmon_index -l A -1 $FQ/$1\_R1.fq.gz -2 $FQ/$1\_R2.fq.gz --validateMappings -o ./salmon_bs/$1 -p 16 --numBootstraps 100