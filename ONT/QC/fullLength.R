source("/wehisan/home/allstaff/d/dong.x/analysis/2020/smchd1/NSC/QC/func.R")
library(ShortRead)
library(GenomicAlignments)
library(dplyr)
library(ggplot2)
library(data.table)
library(GenomicFeatures)
library(Hmisc)
# read and calculate, save to RDS
dir_bam  <- "/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/ONT/transcriptome_mapping"
bams <- list.files(path = dir_bam, pattern="primary.bam$")
for (i in 1:6){
  bam1 <- suppressWarnings(readBam(file.path(dir_bam, bams[[i]])))
  readDF <- makeReadDf(bam1)
  summ <- makeSummaryList(bam1)
  saveRDS(readDF, file = paste("./full_length/", i, ".readDF.RDS"))
  saveRDS(summ, file= paste("./full_length/", i, ".summ.RDS"))
  saveRDS(bam1, file= paste("./full_length/", i, ".bam.RDS"))
  rm(bam1)
  gc()
}

#------------------------
# 
# # prepare for plotting
# 
# # attach tx len from annotation to readDF
# anno <- read.delim("/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/rnasequin_isoforms_2.4.tsv", sep = "\t", stringsAsFactors = FALSE)
# readDF$tx_len <- anno$LENGTH[match(readDF$seqnames, anno$NAME)]
# # remove 'novel' isoforms
# readDF$novelty <- is.na(readDF$tx_len)
# readDF <- readDF[!readDF$novelty,]
# maxLength = max(readDF$tx_len)
# readDF$txLengthGroup <- cut2(readDF$tx_len, cuts = c(0, 500, 1000, 1500,
#                                              2000, maxLength))
# readDF$covFraction <- readDF$width / readDF$tx_len
# library(viridis)
# #calculate some stats by transcript
# txStat <- sapply(unique(readDF$seqnames), function(x){
#   readDF.sel = readDF[readDF$seqnames==x, ]
#   meanCovFrac = mean(readDF.sel$covFraction)
#   medianCovFraction = median(readDF.sel$covFraction)
#   fl95 = sum(readDF.sel$covFraction >= 0.95) / nrow(readDF.sel)
#   fl90 = sum(readDF.sel$covFraction >= 0.90) / nrow(readDF.sel)
#   c(readDF.sel[1, "tx_len"], meanCovFrac, medianCovFraction, fl95, fl90, nrow(readDF.sel))
# }, simplify = TRUE)
# # Fig 3C
# pdf("plots/txLenFL.pdf", height = 5, width = 8)
# ggplot(txStat, aes(x=tx_len, y=fl95, colour = log_count))+
#   scale_x_continuous(trans = "log10") +
#   geom_point() +
#   labs(x = "Annotated transcript length", y = "Fraction of full-length", colour = "log count") +
#   theme_bw() +
#   theme(text = element_text(size = 20)) +
#   scale_colour_viridis()
# dev.off()