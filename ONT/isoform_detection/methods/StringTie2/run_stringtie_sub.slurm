#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --output=stringtie_ont_sub.out
#SBATCH --mail-user=du.m@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

# Author: Mei Du
# Using StringTie2 version 2.1.5 (https://github.com/gpertea/stringtie)

STRINGTIE=/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/stringtie
REF=/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/references/genome_rnasequin_decoychr_2.4_edited_seqremoved.gtf

# running stringtie on each sample
for FILE in *sorted.bam; do
$STRINGTIE/stringtie -L \
-G $REF \
-p 8 \
-A ONT_stringtie_$(basename -- "${FILE%_*}")_sub.txt \
-o ONT_stringtie_$(basename -- "${FILE%_*}")_sub.gtf \
$FILE
done

wait 

# merge multiple stringtie outputs
$STRINGTIE/stringtie --merge \
-G $REF \
-o ONT_stringtie_merged_sub.gtf \
 ONT_stringtie_barcode01_sub.gtf \
 ONT_stringtie_barcode02_sub.gtf \
 ONT_stringtie_barcode03_sub.gtf \
 ONT_stringtie_barcode04_sub.gtf \
 ONT_stringtie_barcode05_sub.gtf \
 ONT_stringtie_barcode06_sub.gtf