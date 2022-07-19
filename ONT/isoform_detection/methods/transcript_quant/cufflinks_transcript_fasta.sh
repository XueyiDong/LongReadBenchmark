# Author: Mei Du
# Using Cufflinks version 2.1.1 (https://github.com/cole-trapnell-lab/cufflinks) 

CUFFLINKS=/stornext/General/data/user_managed/grpu_mritchie_0/Mei/long_read_benchmark/cufflinks
REF=/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/references/genome_rnasequin_decoychr_2.4.fa
ONT_ANNOT=/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/plots/readcount_cor/ONT_annotations

export PATH=$PATH:/stornext/General/data/user_managed/grpu_mritchie_0/Mei/long_read_benchmark/cufflinks

# FLAIR and FLAMES already have fasta

$CUFFLINKS/gffread -w ont_bambu_transcripts.fa -g $REF $ONT_ANNOT/ONT_extended_annotations.gtf

$CUFFLINKS/gffread -w ont_sqanti_transcripts.fa -g $REF $ONT_ANNOT/merged_ont_0.95.collapsed.gff

$CUFFLINKS/gffread -w ont_stringtie_ranscripts.fa -g $REF $ONT_ANNOT/ONT_stringtie_merged_edited.gtf

$CUFFLINKS/gffread -w ont_talon_transcripts.fa -g $REF $ONT_ANNOT/ONT_filtered_annot_talon.gtf


