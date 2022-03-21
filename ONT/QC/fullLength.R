# source("/wehisan/home/allstaff/d/dong.x/analysis/2020/smchd1/NSC/QC/func.R")
suppressPackageStartupMessages({
  library(ShortRead)
  library(GenomicAlignments)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(GenomicFeatures)
  library(Hmisc)
  library(edgeR)
})


#-----------
# # read and calculate, save to RDS
# dir_bam  <- "/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/ONT/transcriptome_mapping/split"
# bams <- list.files(path = dir_bam, pattern=".bam$")
# args = commandArgs(trailingOnly=TRUE)
# i = as.numeric(args[1])
#   bam1 <- suppressWarnings(readBam(file.path(dir_bam, bams[i])))
#   readDF <- makeReadDf(bam1)
#   summ <- makeSummaryList(bam1)
#   saveRDS(readDF, file = paste0("./full_length/test/", bams[i], ".readDF.RDS"))
#   saveRDS(summ, file= paste0("./full_length/test/", bams[i], ".summ.RDS"))
#   saveRDS(bam1, file= paste0("./full_length/test/", bams[i], ".bam.RDS"))
#   rm(bam1)
#   gc()

# #--------------
# # combine readDF for each sample
# dir_rds <- "/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/ONT/QC/full_length"
# args = commandArgs(trailingOnly=TRUE)
# bc = as.numeric(args[1])
# files = list.files(path = dir_rds, pattern="readDF.RDS$")
# smpl = grep(paste0("^barcode0", bc), files)
# rdf = data.frame()
# for(i in smpl){
#   cat(files[i])
#   cat("\n")
#   tmp = readRDS(file.path(dir_rds, files[i]))
#   rdf = rbind(rdf, tmp)
# }
# saveRDS(rdf, file.path(dir_rds, "readDF", paste0("barcode0", bc, ".readDF.RDS")))

#---------------
# # combine all RDF
# dir_rdf <- "/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/ONT/QC/full_length/readDF"
# files = list.files(path = dir_rdf, pattern="readDF.RDS$")
# readDF = data.frame()
# for(i in 1:6){
#     tmp = readRDS(file.path(dir_rdf, files[i]))
#     readDF = rbind(readDF, tmp)
#     rm(tmp)
#     gc()
# }
# saveRDS(rdf, "readDF.RDS")

#------------------------
# # split RDF by tx
# args = commandArgs(trailingOnly=TRUE)
# args = as.numeric(args[1])
# 
# # prepare for plotting
# Sys.time()
# cat("reading RDF.", "\n")
# readDF <- readRDS("readDF2.RDS")
# Sys.time()
# cat("RDF loaded.", "\n")
# 
# tx <- readRDS("tx.RDS")
# tx <- as.character(tx)
# tx.sel <- tx[(1 + 9283 * (args - 1)) : min((9283 * args), length(tx))]
# readDF.sel <- readDF[readDF$seqnames %in% tx.sel, ]
# saveRDS(readDF.sel, paste0("readDF/readDF", args, ".RDS"))

# #--------------
# # calculate tx stats
# # read in single readDF
# dir <- "/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/ONT/QC/readDF"
# args = commandArgs(trailingOnly=TRUE)
# i = as.numeric(args[1])
# RDFs <- list.files(path = dir, pattern="RDS$")
# Sys.time()
# cat("reading RDF", RDFs[i], "\n")
# readDF <- readRDS(file.path(dir, RDFs[i]))
# Sys.time()
# cat("RDF loaded.", "\n")
# gc()
# Sys.time()
# 
# tx <- unique(readDF$seqnames)
# cat("A total of", length(tx),"tx", "\n")
# txStat <- sapply(tx, function(x){
#   # cat("calculating for", substring(x, 1,17), "\n")
#   sel = which(readDF$seqnames==x)
#   meanCovFrac = mean(readDF[sel, ]$covFraction)
#   medianCovFrac = median(readDF[sel, ]$covFraction)
#   fl95 = sum(readDF[sel, ]$covFraction >= 0.95) / length(sel)
#   fl90 = sum(readDF[sel, ]$covFraction >= 0.90) / length(sel)
#   # cat(substring(x, 1,17),  "length:", readDF[sel[1], "tx_len"], ", FL count:", fl95, "\n")
#   return(c(x, readDF[sel[1], "tx_len"], meanCovFrac, medianCovFrac, fl95, fl90, length(sel)))
# }, simplify=TRUE)
# Sys.time()
# cat("saving intermediate output result.", "\n")
# saveRDS(txStat, paste0("txStat/byTx/txStatRaw", i, ".RDS"))
# rm(readDF)
# cat("clean memory.", "\n")
# gc()
# Sys.time()
# cat("transformming tx stats.", "\n")
# txStat <- as.data.frame(t(txStat))
# colnames(txStat) <- c("tx", "tx_len", "mean", "median", "fl95", "fl90", "count")
# gc()
# Sys.time()
# cat("calculation completed.", "\n")
# saveRDS(txStat, paste0("txStat/byTx/txStat", i, ".RDS"))

# -------------------------------------
dir <- "txStat/byTx"
ts <- list.files(dir, "txStat")
txStat <- data.frame()
for(i in 1:17){
  tmp <- readRDS(file.path(dir, ts[[i]]))
  txStat <- rbind(txStat, tmp)
}
rm(tmp)
txStat <- txStat[-1, ]
# ---------------------------------------
# making plots
# gc()
cat("making plot.", "\n")
library(viridis)
pdf("plots/txLenFL.pdf", height = 5, width = 8)
ggplot(txStat, aes(x=tx_len, y=fl95, colour = log_count))+
  scale_x_continuous(trans = "log10") +
  geom_point() +
  labs(x = "Annotated transcript length", y = "Fraction of full-length", colour = "log count") +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  scale_colour_viridis()
ggplot(txStat, aes(x=tx_len, y=fl95))+
  scale_x_continuous(trans = "log10") +
  stat_binhex() +
  theme_bw() +
  labs(x = "Annotated transcript length", y = "Fraction of full-length") +
  scale_fill_viridis(trans = "log10")+
  theme(text=element_text(size = 20))
dev.off()
Sys.time()
