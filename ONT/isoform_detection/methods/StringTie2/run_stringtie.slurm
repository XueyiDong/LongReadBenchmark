#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=stringtie_ont.out
#SBATCH --mail-user=du.m@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

# Author: Mei Du
# Using StringTie2 version 2.1.5 (https://github.com/gpertea/stringtie)

STRINGTIE=/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/stringtie
REF=/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/references/genome_rnasequin_decoychr_2.4_edited.gtf

# running stringtie on each sample
for FILE in *sorted.bam; do
$STRINGTIE/stringtie -L \
-G $REF \
-p 8 \
-o ONT_stringtie_$(basename -- "${FILE%_*}").gtf \
$FILE
done

wait 

# merge multiple stringtie outputs
$STRINGTIE/stringtie --merge \
-G $REF \
-o ONT_stringtie_merged.gtf \
 ONT_stringtie_barcode01.gtf \
 ONT_stringtie_barcode02.gtf \
 ONT_stringtie_barcode03.gtf \
 ONT_stringtie_barcode04.gtf \
 ONT_stringtie_barcode05.gtf \
 ONT_stringtie_barcode06.gtf
