#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --output=salmon_ont_stringtie.out
#SBATCH --mail-user=du.m@wehi.edu.au
#SBATCH --mail-type=END,FAIL

# Author: Mei Du
# Using Salmon version 1.5.2 (https://github.com/COMBINE-lab/salmon)
# mapping-based mode

SALMON=/stornext/General/data/user_managed/grpu_mritchie_0/Mei/long_read_benchmark/salmon
REF=/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/references
ONT_ANNOT=/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/plots/readcount_cor/ONT_annotations
CUFFLINKS=/stornext/General/data/user_managed/grpu_mritchie_0/Mei/long_read_benchmark/cufflinks
FQ=/stornext/Projects/promethion/promethion_access/lab_ritchie/transcr_bench/long_term/transcr_bench/guppy_4_rebasecall_Oct20/barcode01-06

export PATH=$PATH:/stornext/General/data/user_managed/grpu_mritchie_0/Mei/long_read_benchmark/cufflinks
$CUFFLINKS/gffread -w $CUFFLINKS/ont_stringtie_transcripts_new.fa -g $REF/genome_rnasequin_decoychr_2.4.fa $ONT_ANNOT/ONT_stringtie_merged_edited.gtf

# prepare transcriptome index (decoy sequences)
$SALMON/SalmonTools/scripts/generateDecoyTranscriptome.sh \
-b /stornext/System/data/apps/bedtools/bedtools-2.26.0/bin/bedtools \
-m $SALMON/mashmap/mashmap \
-a $ONT_ANNOT/ONT_stringtie_merged_edited.gtf \
-g $REF/genome_rnasequin_decoychr_2.4.fa \
-t $CUFFLINKS/ont_stringtie_transcripts_new.fa \
-o ont_stringtie_new

cd $SALMON/ont_stringtie_new

# transcriptome index
$SALMON/bin/salmon index -t gentrome.fa -i transcripts_index --decoys decoys.txt -k 31

# ONT quantification
$SALMON/bin/salmon quant -i transcripts_index \
-l U \
-r $FQ/barcode01.fq.gz $FQ/barcode02.fq.gz $FQ/barcode03.fq.gz $FQ/barcode04.fq.gz $FQ/barcode05.fq.gz $FQ/barcode06.fq.gz \
--validateMappings \
--threads 16 \
-o transcript_quant

# illumina quantification
ILLUMINA=/stornext/General/data/user_managed/grpu_mritchie_1/long_RNA_benchmark/Illumina_April20/data/raw

$SALMON/bin/salmon quant -i transcripts_index \
-l U \
-r $ILLUMINA/H1975-1_S9_R1.fastq.gz $ILLUMINA/H1975-1_S9_R2.fastq.gz $ILLUMINA/H1975-2_S10_R1.fastq.gz $ILLUMINA/H1975-2_S10_R2.fastq.gz $ILLUMINA/H1975-3_S11_R1.fastq.gz $ILLUMINA/H1975-3_S11_R2.fastq.gz $ILLUMINA/HCC827-1_S12_R1.fastq.gz $ILLUMINA/HCC827-1_S12_R2.fastq.gz $ILLUMINA/HCC827-2_S13_R1.fastq.gz $ILLUMINA/HCC827-2_S13_R2.fastq.gz $ILLUMINA/HCC827-2_S13_topup_R1.fastq.gz $ILLUMINA/HCC827-2_S13_topup_R2.fastq.gz $ILLUMINA/HCC827-5-repeat_S14_R1.fastq.gz $ILLUMINA/HCC827-5-repeat_S14_R2.fastq.gz \
--validateMappings \
--threads 16 \
-o illumina_transcript_quant