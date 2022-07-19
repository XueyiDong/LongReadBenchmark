# bambu pipeline
# Author: Mei Du

BiocManager::install(c("bambu", "Rsamtools", "BSGenome"))
library(Rsamtools)
library(bambu)
library(BSgenome)

# genome
# create fai file first with samtools
gencode_seq = FaFile("/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/references/genome_rnasequin_decoychr_2.4.fa")

# annotations
gencode_gtf = prepareAnnotations("/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/references/genome_rnasequin_decoychr_2.4_edited.gtf")

# BAM files - ONT
bam_loc <- "/stornext/Projects/promethion/promethion_access/lab_ritchie/transcr_bench/long_term/transcr_bench/guppy_4_rebasecall_Oct20_aligned_minimap2"
for (i in (1:6)) {
  i <- paste0(0,i)
  assign(paste0("b",i),BamFile(paste0(bam_loc,"/barcode",i,".bam")))}
merged_ont <- BamFileList(b01, b02, b03, b04, b05, b06)

ont_iso_analysis = bambu(reads = merged_ont, annotations = gencode_gtf, genome = gencode_seq, ncore = 8)

writeBambuOutput(ont_iso_analysis, path = "/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/bambu/outputs", prefix = "ONT_")