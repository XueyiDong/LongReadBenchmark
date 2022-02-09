#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --output=flames_ont_subsample.out
#SBATCH --mail-user=du.m@wehi.edu.au
#SBATCH --mail-type=END,FAIL

# Author: Mei Du
# Using FLAMES version 0.1.0 (https://github.com/LuyiTian/FLAMES)

module load anaconda3
source activate 
conda activate FLAMES

REF=/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/references
FQ_DIR=/stornext/General/data/user_managed/grpu_mritchie_0/Mei/long_read_benchmark/ONT_sub/fastq

/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/FLAMES/python/bulk_long_pipeline.py \
    --gff3 $REF/genome_rnasequin_decoychr_2.4_edited_seqremoved.gtf\
    --fq_dir $FQ_DIR \
    --genomefa $REF/genome_rnasequin_decoychr_2.4.fa \
    --config_file /stornext/General/data/user_managed/grpu_mritchie_0/Mei/long_read_benchmark/FLAMES/config.json \
    --outdir /vast/scratch/users/du.m/FLAMES_sub

