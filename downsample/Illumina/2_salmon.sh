#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --job-name=Illumina_salmon
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=END,FAIL



module load stornext
export PATH=/home/users/allstaff/dong.x/Programs/salmon/salmon-latest_linux_x86_64/bin:$PATH

TX=/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/HumanSequins/gencode.v33.sequins.transcripts.fa
FQ=/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/downsample/Illumina/fastq
# 1 sample 2 number of reads
mkdir -p ./salmon/$1\_$2
salmon quant -i /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/HumanSequins/salmon_index -l A -1 $FQ/$1\_R1.$2.fq.gz -2 $FQ/$1\_R2.$2.fq.gz --validateMappings -o ./salmon/$1\_$2 -p 16 --numBootstraps 100
snrmdiskcopy $FQ/$1\_R1.$2.fq.gz
snrmdiskcopy $FQ/$1\_R2.$2.fq.gz