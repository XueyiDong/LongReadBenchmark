#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=250G
#SBATCH --output=flames_ont_quant.out
#SBATCH --mail-user=du.m@wehi.edu.au
#SBATCH --mail-type=END,FAIL

# Author: Mei Du
# Using FLAMES version 0.1.0 (https://github.com/LuyiTian/FLAMES)
# Run on the full dataset 

module load anaconda3
source activate 
conda activate FLAMES

REF=/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/references
FQ_DIR=/stornext/Projects/promethion/promethion_access/lab_ritchie/transcr_bench/long_term/transcr_bench/guppy_4_rebasecall_Oct20/barcode01-06

# ran in 2 stages due to resource limit
# step 1 - get aligned bam
/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/FLAMES/python/bulk_long_pipeline.py \
    --gff3 $REF/genome_rnasequin_decoychr_2.4_edited.gtf \
    --fq_dir $FQ_DIR
    --genomefa $REF/genome_rnasequin_decoychr_2.4.fa \
    --config_file /stornext/General/data/user_managed/grpu_mritchie_0/Mei/long_read_benchmark/FLAMES/config.json \
    --outdir /vast/scratch/users/du.m/FLAMES


# steps 2-4 using aligned bam as input
/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/FLAMES/python/bulk_long_pipeline_fq_exists.py \
    --gff3 $REF/genome_rnasequin_decoychr_2.4_edited.gtf \
    --inbam /vast/scratch/users/du.m/FLAMES/align2genome.bam \
    --genomefa $REF/genome_rnasequin_decoychr_2.4.fa \
    --config_file /stornext/General/data/user_managed/grpu_mritchie_0/Mei/long_read_benchmark/FLAMES/config_bam.json \
    --outdir /vast/scratch/users/du.m/FLAMES