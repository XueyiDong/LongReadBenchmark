library(ggplot2)
library(Rsubread)


source("src/helper.functions.R")
library(Rsubread)
library(GenomicRanges)
library(plyranges)
library(dplyr)
# read in GTF annotation
anno <- read_gff("../../annotation/HumanSequins/gencode.v33.sequins.gtf")

# Make feature-level annotation
# Get full gene lengths
gene.anno <- anno %>% filter(type == "gene") %>% reduce
# Get intergenic regions as the gaps between the genes
intergenic.anno <- gaps(gene.anno)
# saveRDS(intergenic.anno, "intergenic.anno.RDS")
# Get exonic regions
exons.anno <- anno %>% filter(type == "exon") %>% reduce
# Get the gaps between the exons - these will include introns and intergenic regions
introns.anno <- gaps(exons.anno)
# Check for intergenic regions in the introns object and remove them
summary(introns.anno %in% intergenic.anno)
introns.anno <- introns.anno[!(introns.anno %in% intergenic.anno)]
# Transform the GRanges objects into featureCounts compatible objects
gene.anno <- as.data.frame(gene.anno)
intergenic.anno <- as.data.frame(intergenic.anno)
exons.anno <- as.data.frame(exons.anno)
introns.anno <- as.data.frame(introns.anno)
names(gene.anno) <- names(intergenic.anno) <- names(exons.anno) <- names(introns.anno) <- c("chr", "start", "end", "width", "strand")
gene.anno$GeneID <- rownames(gene.anno)
exons.anno$GeneID <- rownames(exons.anno)
introns.anno$GeneID <- rownames(introns.anno)
intergenic.anno$GeneID <- rownames(intergenic.anno)

### Assign the reads at feature-level 
bam.long <- list.files("/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/ONT/bam/", 
                       pattern = ".bam$", full.names = TRUE)
bam.long <- bam.long[1:6]
bam.short <- list.files("/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/illumina/bam/", 
                        pattern = "sorted.bam$", full.names = TRUE)
bam.short <- bam.short[c(1:5, 8)]
# gene counts
gene.counts.long <- featureCounts(files = bam.long, annot.ext = gene.anno, 
                                  allowMultiOverlap = T,
                                  isLongRead = T,
                                  read2pos = 5,
                                  nthreads = 16,
                                  primaryOnly = TRUE)
saveRDS(gene.counts.long, "counts/gene.counts.long.RDS")
gene.counts.short <- featureCounts(files = bam.short, annot.ext = gene.anno,
                                   allowMultiOverlap = T,
                                   isPairedEnd = TRUE,
                                   read2pos = 5,
                                   nthreads = 16,
                                   primaryOnly = TRUE)
saveRDS(gene.counts.short, "counts/gene.counts.short.RDS")
# exon counts
exon.counts.long <- featureCounts(files = bam.long, annot.ext = exons.anno, 
                                  allowMultiOverlap = T,
                                  isLongRead = T,
                                  read2pos = 5,
                                  nthreads = 16,
                                  primaryOnly = TRUE)
saveRDS(exon.counts.long, "counts/exon.counts.long.RDS")
exon.counts.short <- featureCounts(files = bam.short, annot.ext = exons.anno,
                                   allowMultiOverlap = T,
                                   isPairedEnd = TRUE,
                                   read2pos = 5,
                                   nthreads = 16,
                                   primaryOnly = TRUE)
saveRDS(exon.counts.short, "counts/exon.counts.short.RDS")
# intron counts
intron.counts.long <- featureCounts(files = bam.long, annot.ext = introns.anno, 
                                  allowMultiOverlap = T,
                                  isLongRead = T,
                                  read2pos = 5,
                                  nthreads = 16,
                                  primaryOnly = TRUE)
saveRDS(intron.counts.long, "counts/intron.counts.long.RDS")
intron.counts.short <- featureCounts(files = bam.short, annot.ext = introns.anno,
                                   allowMultiOverlap = T,
                                   isPairedEnd = TRUE,
                                   read2pos = 5,
                                   nthreads = 16,
                                   primaryOnly = TRUE)
saveRDS(intron.counts.short, "counts/intron.counts.short.RDS")
#intergenic counts
intergenic.counts.long <- featureCounts(files = bam.long, annot.ext = intergenic.anno, 
                                  allowMultiOverlap = T,
                                  isLongRead = T,
                                  read2pos = 5,
                                  nthreads = 16,
                                  primaryOnly = TRUE)
saveRDS(intergenic.counts.long, "counts/intergenic.counts.long.RDS")
intergenic.counts.short <- featureCounts(files = bam.short, annot.ext = intergenic.anno,
                                   allowMultiOverlap = T,
                                   isPairedEnd = TRUE,
                                   read2pos = 5,
                                   nthreads = 16,
                                   primaryOnly = TRUE)
saveRDS(intergenic.counts.short, "counts/intergenic.counts.short.RDS")
