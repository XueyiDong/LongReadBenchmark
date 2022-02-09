#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=100G
#SBATCH --job-name=talon
#SBATCH --output=talon_ont_filtering.out
#SBATCH --mail-user=du.m@wehi.edu.au
#SBATCH --mail-type=END,FAIL

# Author: Mei Du
# Using TALON version 5.0.0 (https://github.com/mortazavilab/TALON/)

REF_DIR=/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/references
SAM_DIR=/vast/scratch/users/du.m/TALON

#module load python/3.7.0
module load bedtools
#module load htslib/1.7

# internal priming check
mkdir labeled_ont
for FILE in ${SAM_DIR}/*.sam; do
talon_label_reads --f ${FILE} \
	--g ${REF_DIR}/genome_rnasequin_decoychr_2.4.fa \
	--t 8 \
	--deleteTmp \
	--o labeled_ont/$(basename -- "${FILE%.*}")
done

# database initialisation
talon_initialize_database \
	--f ${REF_DIR}/genome_rnasequin_decoychr_2.4_edited.gtf \
	--a REF_annot \
	--g REF \
	--o ONT_talon

export TMPDIR=$(mktemp -d --tmpdir=/vast/scratch/users/du.m/TALON)
# annotation - run separately, took 4 days, 64 CPU, 200 threads, 600GB
talon \
	--f config.csv \
	--db ONT_talon.db \
	--build REF \
	--o ONT \
	--t 200
rm -rf $TMPDIR
		
# abundance and filtering
talon_summarize \
	--db ONT_talon.db \
	--v \
	--o ONT

# unfiltered
talon_abundance \
	--db ONT_talon.db \
	-a REF_annot \
	--build REF \
	--o ONT_unfiltered

# filtered
talon_filter_transcripts \
	--db ONT_talon.db \
	-a REF_annot \
	--o ONT_filtered_transcripts.csv

talon_abundance \
	--db ONT_talon.db \
	--whitelist ONT_filtered_transcripts.csv \
	-a REF_annot \
	--build REF \
	--o ONT

# GTF
talon_create_GTF \
	--db ONT_talon.db \
	--whitelist ONT_filtered_transcripts.csv \
	-a REF_annot \
	--build REF \
	--o ONT_filtered_annot