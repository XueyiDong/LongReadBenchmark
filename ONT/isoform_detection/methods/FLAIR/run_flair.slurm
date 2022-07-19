#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G
#SBATCH --output=flair_ont.out
#SBATCH --mail-user=du.m@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

# Author: Mei Du
# Using FLAIR version 1.5.0 (https://github.com/BrooksLabUCSC/flair)
# Run on the full dataset

module load anaconda3
source activate
conda activate flair_env


FLAIR_DIR=/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/FLAIR
REF=/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/references
FQ_DIR=/stornext/Projects/promethion/promethion_access/lab_ritchie/transcr_bench/long_term/transcr_bench/guppy_4_rebasecall_Oct20/barcode01-06
SUB_FQ_DIR=/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/ONT_aligned/subsample_fq
OUT=/vast/scratch/users/du.m/FLAIR

# pipeline divided into 2 stages due to resource limit - 1. align, correct 2. collapse, quantify
# align
python $FLAIR_DIR/flair.py align \
-g $REF/genome_rnasequin_decoychr_2.4.fa \
-r $FQ_DIR/barcode01.fq.gz $FQ_DIR/barcode02.fq.gz $FQ_DIR/barcode03.fq.gz $FQ_DIR/barcode04.fq.gz $FQ_DIR/barcode05.fq.gz $FQ_DIR/barcode06.fq.gz \
-t 16 \
-v1.3

# correct 
python $FLAIR_DIR/flair.py correct \
-q $OUT/flair.aligned.bed \
-f $REF/genome_rnasequin_decoychr_2.4_edited.gtf \
-g $REF/genome_rnasequin_decoychr_2.4.fa \
--nvrna \
-t 16

# collapse
python $FLAIR_DIR/flair.py collapse -g $REF/genome_rnasequin_decoychr_2.4.fa \
-f $REF/genome_rnasequin_decoychr_2.4_edited.gtf \
-r $FQ_DIR/barcode01.fq.gz $FQ_DIR/barcode02.fq.gz $FQ_DIR/barcode03.fq.gz $FQ_DIR/barcode04.fq.gz $FQ_DIR/barcode05.fq.gz $FQ_DIR/barcode06.fq.gz \
-q $OUT/flair_all_corrected.bed \
-s 10 \
-t 16 \
--trust_ends \
--temp_dir $OUT/tmp

# quantify
python $FLAIR_DIR/flair.py quantify \
-r $OUT/reads_manifest.tsv \
-i $OUT/flair.collapse.isoforms.fa \
--trust_ends \
-t 16 \
--temp_dir $OUT/tmp
