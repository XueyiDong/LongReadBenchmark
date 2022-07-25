#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

cd /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/illumina

export PATH=/home/users/allstaff/dong.x/Programs/salmon/salmon-latest_linux_x86_64/bin:$PATH

TX=/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/HumanSequins/gencode.v33.sequins.transcripts.fa
FQ=/stornext/General/data/user_managed/grpu_mritchie_1/long_RNA_benchmark/Illumina_April20/data/raw/
for SAMPLE in H1975-1_S9 H1975-2_S10 H1975-3_S11 HCC827-1_S12 HCC827-2_S13 HCC827-5-repeat_S14
do 
# mkdir -p ./salmon/$SAMPLE
# salmon quant -i /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/HumanSequins/salmon_index -l A -1 $FQ/$SAMPLE\_R1.fastq.gz -2 $FQ/$SAMPLE\_R2.fastq.gz --validateMappings -o ./salmon/$SAMPLE
mkdir -p ./salmon_bs/$SAMPLE
salmon quant -i /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/HumanSequins/salmon_index -l A -1 $FQ/$SAMPLE\_R1.fastq.gz -2 $FQ/$SAMPLE\_R2.fastq.gz --validateMappings -o ./salmon_bs/$SAMPLE -p 16 --numBootstraps 100
done