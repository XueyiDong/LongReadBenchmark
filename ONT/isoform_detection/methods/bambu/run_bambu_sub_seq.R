# bambu pipeline using downsampled dataset for sequin analysis
# Author: Mei Du

#BiocManager::install(c("bambu", "Rsamtools", "BSGenome"))
library(Rsamtools)
library(bambu)
library(BSgenome)

# genome
# create fai file first with samtools
gencode_seq = FaFile("/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/references/genome_rnasequin_decoychr_2.4.fa")

# annotations
gencode_gtf = prepareAnnotations("/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/references/genome_rnasequin_decoychr_2.4_edited_seqremoved.gtf")

# BAM files - ONT
sub_bam_loc <- "/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/ONT_aligned/subsample_aligned"
for (i in (1:6)) {
  i <- paste0(0,i)
  assign(paste0("b",i),BamFile(paste0(sub_bam_loc,"/barcode",i,"_sub_10.bam")))}
merged_ont <- BamFileList(b01, b02, b03, b04, b05, b06)

ont_iso_analysis = bambu(reads = merged_ont, annotations = gencode_gtf, genome = gencode_seq, ncore = 16)

writeBambuOutput(ont_iso_analysis, path = "/stornext/General/data/user_managed/grpu_mritchie_0/Mei/long_read_benchmark/ONT_sub/bambu", prefix = "ONT_sub_")