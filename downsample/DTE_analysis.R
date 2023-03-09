# Author: Xueyi Dong
# Aim: explore the performance of DTE analysis under different library size

library(edgeR)
library(limma)
# Read in ONT counts
counts_ONT <- readRDS("./ONT/counts.RDS")
counts_ONT_full <- readRDS("../ONT/counts_ONT_full.RDS")
names(counts_ONT_full) <- paste0(names(counts_ONT_full), ".full")
# check whether row names match
all(rownames(counts_ONT) == rownames(counts_ONT_full))
# doesn't match. need to match row names before merging
m <- match(rownames(counts_ONT_full), rownames(counts_ONT))
counts_ONT <- cbind(counts_ONT[m, ], counts_ONT_full)
# Read in Illumina counts
counts_Illumina <- readRDS("./Illumina/counts.RDS")
counts_Illumina_full <- readRDS("../illumina/counts_Illumina_full.RDS")
names(counts_Illumina_full) <- paste0(names(counts_Illumina_full), "_full")
# check whether row names match
all(rownames(counts_Illumina) == rownames(counts_Illumina_full))
# All match. combine counts from full and downsampled samples.
counts_Illumina <- cbind(counts_Illumina, counts_Illumina_full)
# Keep pure samples only
counts_ONT_pure <- counts_ONT[, grep("barcode|H", colnames(counts_ONT))]
counts_Illumina_pure <- counts_Illumina[, grep("H", colnames(counts_Illumina))]
# Organize ONT sample information
samples_ONT<- as.data.frame(strsplit2(colnames(counts_ONT_pure), "\\."))
colnames(samples_ONT) <- c("sample", "libsize")
samples_ONT$group[samples_ONT$sample %in% paste0("barcode0", 1:3)] <- "H1975"
samples_ONT$group[samples_ONT$sample %in% paste0("barcode0", 4:6)] <- "HCC827"
samples_ONT$group[grep("H", colnames(counts_ONT_pure))] <- strsplit2(colnames(counts_ONT_pure)[grep("H", colnames(counts_ONT_pure))], "\\-")[,1]
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
pdf("plots/MDS_human_all.pdf", height = 4, width = 15)
par(mfrow = c(2, 7))
for(i in 1:length(x_human)){
  cpm.human <- cpm(x_human[[i]], log = TRUE)
  plotMDS(cpm.human, main = names(x_human)[[i]],
          pch = 1,
          col = rep(c("red", "blue"), each = 3))
}
dev.off()
# make MDS plots for sequin transcripts
pdf("plots/MDS_sequin_all.pdf", height = 4, width = 15)
par(mfrow = c(2, 7))
for(i in 1:length(x_sequin)){
  cpm.sequin <- cpm(x_sequin[[i]], log = TRUE)
  plotMDS(cpm.sequin, main = names(x_sequin)[[i]],
          pch = 1,
          col = rep(c("red", "blue"),each = 3))
}
dev.off()
# DE analysis. We keep the topTags output.
tt_human <- lapply(x_human, function(x){
  design <- model.matrix(~x$samples$group)
  x <- estimateDisp(x, design)
  qlfit <- glmQLFit(x, design)
  res <- glmQLFTest(qlfit)
  return(as.data.frame(topTags(res, n = Inf)))
})
tt_sequin<- lapply(x_sequin, function(x){
  design <- model.matrix(~x$samples$group)
  x <- estimateDisp(x, design)
  qlfit <- glmQLFit(x, design)
  res <- glmQLFTest(qlfit)
  return(as.data.frame(topTags(res, n = Inf)))
})
# extract DE transcripts
DE_human <- lapply(tt_human, function(x){
  rownames(x)[x$FDR < 0.05]
})
DE_sequin <- lapply(tt_sequin, function(x){
  rownames(x)[x$FDR < 0.05]
})
