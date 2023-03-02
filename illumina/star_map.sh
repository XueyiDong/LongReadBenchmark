#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

module load STAR/2.7.9a
FADIR=/stornext/General/data/user_managed/grpu_mritchie_1/long_RNA_benchmark/Illumina_April20/data/raw/

mkdir -p ./STAR/$1

STAR --runThreadN 16\
 --genomeDir /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/HumanSequins/GenomeDir\
 --readFilesIn $FADIR/$1\_R1.fastq.gz $FADIR/$1\_R2.fastq.gz\
 --readFilesCommand zcat\
 --outFilterType BySJout\
 --outSAMtype BAM Unsorted\
 --quantMode TranscriptomeSAM\
 --outFileNamePrefix ./STAR/$1/\
 --outSAMattributes NH HI AS NM MD