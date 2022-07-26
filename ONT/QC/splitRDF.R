suppressPackageStartupMessages({
  library(ShortRead)
  library(GenomicAlignments)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(GenomicFeatures)
  library(Hmisc)
})

args = commandArgs(trailingOnly=TRUE)
args = as.numeric(args[1])

# prepare for plotting
Sys.time()
cat("reading RDF.", "\n")
readDF <- readRDS("readDF2.RDS")
Sys.time()
cat("RDF loaded.", "\n")
# split read DF by transcript to avoid too large memory usage
tx <- readRDS("tx.RDS")
tx <- as.character(tx)
tx.sel <- tx[(1 + 9283 * (args - 1)) : min((9283 * args), length(tx))]
readDF.sel <- readDF[readDF$seqnames %in% tx.sel, ]
saveRDS(readDF.sel, paste0("readDF/readDF", args, ".RDS"))