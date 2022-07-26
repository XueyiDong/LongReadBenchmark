#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

cd /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/ONT

export PATH=/home/users/allstaff/dong.x/Programs/salmon/salmon-latest_linux_x86_64/bin:$PATH

module load samtools
TX=/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/HumanSequins/gencode.v33.sequins.transcripts.fa
for SAMPLE in barcode{01..06}
do 
mkdir -p ./salmon_bs/$SAMPLE
samtools view -F 256 -b -@ 8 ./transcriptome_mapping/$SAMPLE.bam > ./transcriptome_mapping/$SAMPLE.primary.bam
salmon quant -t $TX -l A -a ./transcriptome_mapping/$SAMPLE.primary.bam -o ./salmon_bs/$SAMPLE --numBootstraps 100 -p 8
done