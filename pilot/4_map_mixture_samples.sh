#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --job-name=map_pilot
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

cd /stornext/General/data/user_managed/grpu_mritchie_1/Mixture/MainStudy/in_silico_mix
module load subread/1.5.2

# 1 sample name
subread-align -t 0 -T 12 -i /wehisan/home/allstaff/d/dong.x/personal/index/hg19/hg19_index -r $1.fastq.gz -o bam/$1.bam