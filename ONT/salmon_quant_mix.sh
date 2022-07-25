#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

cd /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/ONT/
module load samtools

BAMDIR=/stornext/General/data/user_managed/grpu_mritchie_0/XueyiDong/long_read_benchmark/ONT/mix_fq/bam_tx

export PATH=/home/users/allstaff/dong.x/Programs/salmon/salmon-latest_linux_x86_64/bin:$PATH

TX=/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/HumanSequins/gencode.v33.sequins.transcripts.fa
for SAMPLE in mix_025075_{1..3} mix_050050_{1..3} mix_075025_{1..3}
# SAMPLE=mix_075025_2
do 
mkdir -p ./salmon_bs/$SAMPLE
samtools view -F 256 -b -@ 8 $BAMDIR/$SAMPLE.bam > $BAMDIR/$SAMPLE.primary.bam
salmon quant -t $TX -l A -a $BAMDIR/$SAMPLE.primary.bam -o ./salmon_bs/$SAMPLE --numBootstraps 100 -p 16
done