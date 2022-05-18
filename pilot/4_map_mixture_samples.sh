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

# # rep1
# subread-align -t 0 -T 12 -i /wehisan/home/allstaff/d/dong.x/personal/index/hg19/hg19_index -r /stornext/General/data/user_managed/grpu_mritchie_1/Mixture/MainStudy/R1-025sub.fastq -o /stornext/General/data/user_managed/grpu_mritchie_1/Mixture/MainStudy/R1-025sub.bam
# subread-align -t 0 -T 12 -i /wehisan/home/allstaff/d/dong.x/personal/index/hg19/hg19_index -r /stornext/General/data/user_managed/grpu_mritchie_1/Mixture/MainStudy/R1-050sub.fastq -o /stornext/General/data/user_managed/grpu_mritchie_1/Mixture/MainStudy/R1-050sub.bam
# subread-align -t 0 -T 12 -i /wehisan/home/allstaff/d/dong.x/personal/index/hg19/hg19_index -r /stornext/General/data/user_managed/grpu_mritchie_1/Mixture/MainStudy/R1-075sub.fastq -o /stornext/General/data/user_managed/grpu_mritchie_1/Mixture/MainStudy/R1-075sub.bam
# # rep2
# subread-align -t 0 -T 12 -i /wehisan/home/allstaff/d/dong.x/personal/index/hg19/hg19_index -r /stornext/General/data/user_managed/grpu_mritchie_1/Mixture/MainStudy/R2-025sub.fastq -o /stornext/General/data/user_managed/grpu_mritchie_1/Mixture/MainStudy/R2-025sub.bam
# subread-align -t 0 -T 12 -i /wehisan/home/allstaff/d/dong.x/personal/index/hg19/hg19_index -r /stornext/General/data/user_managed/grpu_mritchie_1/Mixture/MainStudy/R2-050sub.fastq -o /stornext/General/data/user_managed/grpu_mritchie_1/Mixture/MainStudy/R2-050sub.bam
# subread-align -t 0 -T 12 -i /wehisan/home/allstaff/d/dong.x/personal/index/hg19/hg19_index -r /stornext/General/data/user_managed/grpu_mritchie_1/Mixture/MainStudy/R2-075sub.fastq -o /stornext/General/data/user_managed/grpu_mritchie_1/Mixture/MainStudy/R2-075sub.bam
# # rep3
# subread-align -t 0 -T 12 -i /wehisan/home/allstaff/d/dong.x/personal/index/hg19/hg19_index -r /stornext/General/data/user_managed/grpu_mritchie_1/Mixture/MainStudy/R3-025sub.fastq -o /stornext/General/data/user_managed/grpu_mritchie_1/Mixture/MainStudy/R3-025sub.bam
# subread-align -t 0 -T 12 -i /wehisan/home/allstaff/d/dong.x/personal/index/hg19/hg19_index -r /stornext/General/data/user_managed/grpu_mritchie_1/Mixture/MainStudy/R3-050sub.fastq -o /stornext/General/data/user_managed/grpu_mritchie_1/Mixture/MainStudy/R3-050sub.bam
# subread-align -t 0 -T 12 -i /wehisan/home/allstaff/d/dong.x/personal/index/hg19/hg19_index -r /stornext/General/data/user_managed/grpu_mritchie_1/Mixture/MainStudy/R3-075sub.fastq -o /stornext/General/data/user_managed/grpu_mritchie_1/Mixture/MainStudy/R3-075sub.bam