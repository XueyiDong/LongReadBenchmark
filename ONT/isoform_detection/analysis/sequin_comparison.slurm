#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=sequin_gffcompare.out
#SBATCH --mail-user=du.m@wehi.edu.au
#SBATCH --mail-type=END,FAIL

GFFCOMPARE=/stornext/General/data/user_managed/grpu_mritchie_0/Mei/long_read_benchmark/gffcompare/gffcompare
ANNOT=/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/plots/readcount_cor/ONT_annotations
SUB_ANNOT=/stornext/General/data/user_managed/grpu_mritchie_0/Mei/long_read_benchmark/ONT_sub/downsampled_analysis
REF=/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/references/genome_rnasequin_decoychr_2.4_edited.gtf 

# bambu
$GFFCOMPARE \
$ANNOT/ONT_extended_annotations_edited.gtf \
$SUB_ANNOT/bambu/ONT_sub_extended_annotations.gtf \
-o bambu \
-r $REF \
-V \
--debug

# FLAIR
$GFFCOMPARE \
$ANNOT/flair.collapse.isoforms.gtf \
$SUB_ANNOT/FLAIR/flair.collapse.isoforms.gtf \
-o flair \
-r $REF \
-V \
--debug

# FLAMES
$GFFCOMPARE \
$ANNOT/isoform_annotated.filtered_flames.gff3 \
$SUB_ANNOT/FLAMES/isoform_annotated.filtered.gff3 \
-o flames \
-r $REF \
-V \
--debug

# SQANTI3
$GFFCOMPARE \
$ANNOT/merged_ont.collapsed_corrected.gtf \
$SUB_ANNOT/SQANTI3/merged_ont_sub.collapsed_corrected.gtf \
-o sqanti \
-r $REF \
-V \
--debug

# StringTie
$GFFCOMPARE \
$ANNOT/ONT_stringtie_merged_edited.gtf \
$SUB_ANNOT/stringtie/ONT_stringtie_merged_sub.gtf \
-o stringtie \
-r $REF \
-V \
--debug

# TALON
$GFFCOMPARE \
$ANNOT/ONT_filtered_annot_talon.gtf \
$SUB_ANNOT/TALON/ONT_sub_filtered_annot_talon.gtf \
-o talon \
-r $REF \
-V \
--debug