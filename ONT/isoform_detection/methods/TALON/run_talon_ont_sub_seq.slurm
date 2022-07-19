#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --job-name=talon
#SBATCH --output=talon_ont_sub.out
#SBATCH --mail-user=du.m@wehi.edu.au
#SBATCH --mail-type=END,FAIL

# Author: Mei Du
# Using TALON version 5.0.0 (https://github.com/mortazavilab/TALON/)
# Run on the downsampled dataset with the modified GTF

REF_DIR=/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/references
SAM_DIR=/vast/scratch/users/du.m/subsample/TALON

#module load python/3.7.0
module load bedtools
#module load htslib/1.7

# database initialisation
talon_initialize_database \
	--f ${REF_DIR}/genome_rnasequin_decoychr_2.4_edited_seqremoved.gtf \
	--a REF_annot \
	--g REF \
	--o ONT_sub_talon

mkdir -p labeled_ont

# internal priming check
for FILE in *clean.sam; do
talon_label_reads --f ${FILE} \
	--g ${REF_DIR}/genome_rnasequin_decoychr_2.4.fa \
	--t 8 \
	--ar 20 \
	--deleteTmp \
	--o labeled_ont/$(basename -- "${FILE%.*}")
done

export TMPDIR=$(mktemp -d --tmpdir=/vast/scratch/users/du.m/subsample/TALON)
# annotation
talon \
	--f config_ont_sub.csv \
	--db ONT_sub_talon.db \
	--build REF \
	--o ONT_sub \
	--t 8
rm -rf $TMPDIR

# abundance and filtering
talon_summarize \
	--db ONT_sub_talon.db \
	--v \
	--o ONT_sub

# unfiltered
talon_abundance \
	--db ONT_sub_talon.db \
	-a REF_annot \
	--build REF \
	--o ONT_sub_unfiltered

# filtered
talon_filter_transcripts \
	--db ONT_sub_talon.db \
	-a REF_annot \
	--o ONT_sub_filtered_transcripts.csv

talon_abundance \
	--db ONT_sub_talon.db \
	--whitelist ONT_sub_filtered_transcripts.csv \
	-a REF_annot \
	--build REF \
	--o ONT_sub

# GTF
talon_create_GTF \
	--db ONT_sub_talon.db \
	--whitelist ONT_sub_filtered_transcripts.csv \
	-a REF_annot \
	--build REF \
	--o ONT_sub_filtered_annot