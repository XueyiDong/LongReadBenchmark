# Author: Xueyi Dong
# Aim: explore the performance of DTE analysis under different library size

library(edgeR)
library(limma)
# Read in counts and keep only pure samples
counts_ONT <- readRDS("./ONT/counts.RDS")
counts_Illumina <- readRDS("./Illumina/counts.RDS")
counts_ONT_pure <- counts_ONT[, grep("barcode", colnames(counts_ONT))]
counts_Illumina_pure <- counts_Illumina[, grep("H", colnames(counts_Illumina))]
# Organize ONT sample information
samples_ONT<- as.data.frame(strsplit2(colnames(counts_ONT_pure), "\\."))
colnames(samples_ONT) <- c("sample", "libsize")
samples_ONT$group[samples_ONT$sample %in% paste0("barcode0", 1:3)] <- "H1975"
samples_ONT$group[samples_ONT$sample %in% paste0("barcode0", 4:6)] <- "HCC827"
rownames(samples_ONT) <- colnames(counts_ONT_pure)
# Organize Illumina sample information
samples_Illumina <- as.data.frame(strsplit2(colnames(counts_Illumina_pure), "_"))
samples_Illumina[!(samples_Illumina[, 4] == ""), 3] <- samples_Illumina[!(samples_Illumina[, 4] == ""), 4]
samples_Illumina <- samples_Illumina[, c(1, 3)]
colnames(samples_Illumina) <- c("sample", "libsize")
samples_Illumina$group <- strsplit2(samples_Illumina$sample, "-")[,1]
rownames(samples_Illumina) <- colnames(counts_Illumina_pure)
# For each library size: make a DGEList object, preprocess data, separate human and sequins transcripts, do DE analysis H1975 vs HCC827
# make ONT DGEList and filter out lowly expressed transcripts for each library size
x_ONT <- lapply(unique(samples_ONT$libsize), function(x){
  ob <- DGEList(counts = counts_ONT_pure[, samples_ONT$libsize == x],
                samples = samples_ONT[samples_ONT$libsize == x, ])
  filter <- filterByExpr(ob)
  ob <- ob[filter, , keep.lib.sizes = FALSE]
  return(ob)
})
names(x_ONT) <- paste0("ONT_", unique(samples_ONT$libsize))
# make Illumina DGEList and filter out lowly expressed transcripts for each library size
x_Illumina <- lapply(unique(samples_Illumina$libsize), function(x){
  ob <- DGEList(counts = counts_Illumina_pure[, samples_Illumina$libsize == x],
                samples = samples_Illumina[samples_Illumina$libsize == x, ])
  filter <- filterByExpr(ob)
  ob <- ob[filter, , keep.lib.sizes = FALSE]
  return(ob)
})
names(x_Illumina) <- paste0("Illumina_", unique(samples_Illumina$libsize))
# combine lists
x <- append(x_ONT, x_Illumina)
# separate human and sequin transcripts and normalize
x_human <- lapply(x, function(x){
  human <- grepl("^ENST", rownames(x))
  x_human <- x[human, ]
  return(calcNormFactors(x_human))
})
x_sequin <- lapply(x, function(x){
  sequin <- grepl("^R", rownames(x))
  x_sequin <- x[sequin, ]
  return(calcNormFactors(x_sequin))
})
# make MDS plots for human transcripts

