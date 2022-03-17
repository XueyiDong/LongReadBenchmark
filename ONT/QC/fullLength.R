# source("/wehisan/home/allstaff/d/dong.x/analysis/2020/smchd1/NSC/QC/func.R")
library(ShortRead)
library(GenomicAlignments)
library(dplyr)
library(ggplot2)
library(data.table)
library(GenomicFeatures)
library(Hmisc)
# library(edgeR)

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

args = commandArgs(trailingOnly=TRUE)
args = as.numeric(args[1])

# prepare for plotting
Sys.time()
cat("reading RDF.", "\n")
readDF <- readRDS("readDF2.RDS")
Sys.time()
cat("RDF loaded.", "\n")
# gc()
# readDF$seqnames <- as.character(readDF$seqnames)
# readDF$rname <- as.character(readDF$rname)
# # get tx len and attach to readDF
# cat("reading DGE", "\n")
# dge <- readRDS("../../longvsshort/dge.rds")
# gc()
# readDF$tx_len <- dge$genes$Length[match(readDF$seqnames, rownames(dge))]
# # remove 'novel' isoforms
# # readDF$novelty <- is.na(readDF$tx_len)
# # readDF <- readDF[!readDF$novelty,]
# maxLength = max(readDF$tx_len)
# readDF$txLengthGroup <- cut2(readDF$tx_len, cuts = c(0, 500, 1000, 1500,
#                                              2000, maxLength))
# readDF$covFraction <- readDF$width / readDF$tx_len
# cat("Saving RDF.", "\n")
# saveRDS(readDF, "readDF2.RDS")
# gc()
library(parallel)
#calculate some stats by transcript
Sys.time()
cat("calculating tx stats for tx", 1 + 9283 * (args - 1), "to", 9283 * args, "\n")
nCores = detectCores()
cat("number of cores:", nCores, "\n")
# tx <- unique(readDF$seqnames)
# cat("number of transcripts:", length(tx), "\n")
# saveRDS(tx, "tx.RDS")
tx <- readRDS("tx.RDS")
tx <- as.character(tx)
txStat <- mclapply(tx[(1 + 9283 * (args - 1)) : 9283 * args], function(x){
  cat("calculating for tx", x, "\n")
  readDF.sel = readDF[readDF$seqnames==x, ]
  meanCovFrac = mean(readDF.sel$covFraction)
  medianCovFrac = median(readDF.sel$covFraction)
  fl95 = sum(readDF.sel$covFraction >= 0.95) / nrow(readDF.sel)
  fl90 = sum(readDF.sel$covFraction >= 0.90) / nrow(readDF.sel)
  cat("Tx length:", readDF.sel[1, "tx_len"], "full length fraction:", fl95, "\n")
  len.tmp = readDF.sel[1, "tx_len"]
  rm(readDF.sel)
  return(c(len.tmp, meanCovFrac, medianCovFrac, fl95, fl90, nrow(readDF.sel)))
}, mc.cores = nCores)
Sys.time()
cat("saving intermediate output result.", "\n")
saveRDS(txStat, paste0("txStat", args, ".RDS"))
rm(readDF)
cat("clean memory.", "\n")
gc()
Sys.time()
cat("transformming tx stats.", "\n")
ind <- sapply(txStat, is.null, simplify=TRUE)
table(ind)
txStat <- matrix(unlist(txStat), nrow = 6)
txStat <- as.data.frame(t(txStat))
colnames(txStat) <- c("tx_len", "mean", "median", "fl95", "fl90", "count")
txStat$log_count <- log(txStat$count)
txStat$tx <- tx[(1 + 9283 * (args - 1)) : 9283 * args][which(ind==FALSE)]
dge <- readRDS("../../longvsshort/dge.rds")
dge$genes$totalCount <- rowSums(dge$counts[,1:6])
txStat$count_dge <- dge$genes$totalCount[match(txStat$tx, rownames(dge))]
# txStat$transcript <- tx
gc()
Sys.time()
cat("calculation completed.", "\n")
saveRDS(txStat, paste0("txStat_final", args, ".RDS"))
# gc()
# cat("making plot.", "\n")
# Fig 3C
# library(viridis)
# pdf("plots/txLenFL.pdf", height = 5, width = 8)
# ggplot(txStat, aes(x=tx_len, y=fl95, colour = log_count))+
#   scale_x_continuous(trans = "log10") +
#   geom_point() +
#   labs(x = "Annotated transcript length", y = "Fraction of full-length", colour = "log count") +
#   theme_bw() +
#   theme(text = element_text(size = 20)) +
#   scale_colour_viridis()
# ggplot(txStat, aes(x=tx_len, y=fl95))+
#   scale_x_continuous(trans = "log10") +
#   stat_binhex() +
#   theme_bw() +
#   labs(x = "Annotated transcript length", y = "Fraction of full-length") +
#   scale_fill_viridis(trans = "log10")+
#   theme(text=element_text(size = 20)) 
# dev.off()
Sys.time()
