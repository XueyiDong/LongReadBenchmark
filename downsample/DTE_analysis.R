# Author: Xueyi Dong
# Aim: explore the performance of DTE analysis under different library size

library(edgeR)
library(limma)
library(ggplot2)
library(UpSetR)
# Read in ONT counts
counts_ONT <- readRDS("./ONT/counts.RDS")
counts_ONT_full <- readRDS("../ONT/counts_ONT_full.RDS")
names(counts_ONT_full) <- paste0(names(counts_ONT_full), ".full")
# check whether row names match
all(rownames(counts_ONT) == rownames(counts_ONT_full))
# All match. combine counts from full and downsampled samples
counts_ONT <- cbind(counts_ONT, counts_ONT_full)
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
# samples_ONT$group[grep("H", colnames(counts_ONT_pure))] <- strsplit2(colnames(counts_ONT_pure)[grep("H", colnames(counts_ONT_pure))], "\\-")[,1]
# rownames(samples_ONT) <- colnames(counts_ONT_pure)
samples_ONT$libsize[nchar(samples_ONT$libsize) == 2] <- paste0("0", samples_ONT$libsize[nchar(samples_ONT$libsize) == 2])
# Organize Illumina sample information
samples_Illumina <- as.data.frame(strsplit2(colnames(counts_Illumina_pure), "_"))
samples_Illumina[!(samples_Illumina[, 4] == ""), 3] <- samples_Illumina[!(samples_Illumina[, 4] == ""), 4]
samples_Illumina <- samples_Illumina[, c(1, 3)]
colnames(samples_Illumina) <- c("sample", "libsize")
samples_Illumina$group <- strsplit2(samples_Illumina$sample, "-")[,1]
rownames(samples_Illumina) <- colnames(counts_Illumina_pure)
samples_Illumina$libsize[nchar(samples_Illumina$libsize) == 2] <- paste0("0", samples_Illumina$libsize[nchar(samples_Illumina$libsize) == 2])
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
x <- x[order(names(x))]
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
  plotMDS(cpm.human, main = names(x_human)[i],
          xlim = c(-3.4, 3.4),
          ylim = c(-3.6, 4.1),
          pch = 1,
          col = rep(c("red", "blue"), each = 3))
}
dev.off()
# make MDS plots for sequin transcripts
pdf("plots/MDS_sequin_all.pdf", height = 4, width = 15)
par(mfrow = c(2, 7))
for(i in 1:length(x_sequin)){
  cpm.sequin <- cpm(x_sequin[[i]], log = TRUE)
  plotMDS(cpm.sequin, main = names(x_sequin)[i],
          xlim = c(-1.4, 1.4),
          ylim = c(-0.36, 0.39),
          pch = 1,
          col = rep(c("red", "blue"),each = 3))
}
dev.off()
# DE analysis using edgeR quasiLikelihood
qlres_human <- lapply(x_human, function(x){
  design <- model.matrix(~x$samples$group)
  x <- estimateDisp(x, design)
  qlfit <- glmQLFit(x, design)
  res <- glmQLFTest(qlfit)
  return(res)
})
qlres_sequin <- lapply(x_sequin, function(x){
  design <- model.matrix(~x$samples$group)
  x <- estimateDisp(x, design)
  qlfit <- glmQLFit(x, design)
  res <- glmQLFTest(qlfit)
  return(res)
})
# MD plots
pdf("plots/MD_human.pdf", height = 4, width = 15)
par(mfrow = c(2, 7))
for(i in 1:length(qlres_human)){
  plotMD(qlres_human[[i]], status = decideTestsDGE(qlres_human[[i]]),
         values = c(1, -1), col = c("red", "blue"),
         main = names(qlres_human)[i])
}
dev.off()
pdf("plots/MD_sequin.pdf", height = 4, width = 15)
par(mfrow = c(2, 7))
for(i in 1:length(qlres_sequin)){
  plotMD(qlres_sequin[[i]], status = decideTestsDGE(qlres_sequin[[i]]),
         values = c(1, -1), col = c("red", "blue"),
         main = names(qlres_sequin)[i])
}
dev.off()
# topTags output
tt_human <- lapply(qlres_human, function(x){
  return(as.data.frame(topTags(x, n = Inf)))
})
tt_sequin<- lapply(qlres_sequin, function(x){
  return(as.data.frame(topTags(x, n = Inf)))
})
# extract DE transcripts
DE_human <- lapply(tt_human, function(x){
  rownames(x)[x$FDR < 0.05]
})
DE_sequin <- lapply(tt_sequin, function(x){
  rownames(x)[x$FDR < 0.05]
})
# Plot power for human
power_human <- sapply(DE_human, length, simplify = TRUE)
res_human <- data.frame(
  power = power_human,
  dataset = strsplit2(names(power_human), "_")[, 1],
  libsize = strsplit2(names(power_human), "_")[, 2]
)
res_human$libsize = factor(res_human$libsize, levels = c(paste0(c("01", "03", "05", "10", "15", "20"), "M"), "full"))
pdf("plots/power_human.pdf", height = 5, width = 8)
ggplot(res_human, aes(x = libsize, y = power, colour = dataset, group = dataset)) +
  geom_line() +
  geom_point() +
  scale_colour_manual(values = c("#FCB344", "#438DAC")) +
  labs(x = "Number of reads", y = "Number of DE transcripts", colour = "Dataset") +
  theme_bw() +
  theme(text = element_text(size = 20))
dev.off()
# Plot power for sequins
power_sequin <- sapply(DE_sequin, length, simplify = TRUE)
res_sequin <- data.frame(
  power = power_sequin,
  dataset = strsplit2(names(power_sequin), "_")[, 1],
  libsize = strsplit2(names(power_sequin), "_")[, 2]
)
res_sequin$libsize <- factor(res_sequin$libsize, levels = c(paste0(c("01", "03", "05", "10", "15", "20"), "M"), "full"))
pdf("plots/power_sequin.pdf", height = 5, width = 8)
ggplot(res_sequin, aes(x = libsize, y = power, colour = dataset, group = dataset)) +
  geom_line() +
  geom_point() +
  scale_colour_manual(values = c("#FCB344", "#438DAC")) +
  labs(x = "Number of reads", y = "Number of DE transcripts", colour = "Dataset") +
  theme_bw() +
  theme(text = element_text(size = 20))
dev.off()
# upset plot
pdf("plots/upset_human_ONT.pdf", height = 5, width = 12)
upset(fromList(DE_human[1:7]), nsets = 7, order.by = "freq",
      text.scale = c(1.5, 1.5, 1.5, 1.2, 1.2, 1.5))
dev.off()
pdf("plots/upset_human_Illumina.pdf", height = 5, width = 12)
upset(fromList(DE_human[8:14]), nsets = 7, order.by = "freq",
      text.scale = c(1.5, 1.5, 1.5, 1.2, 1.2, 1.5))
dev.off()
pdf("plots/upset_human.pdf", height = 7, width = 15)
upset(fromList(DE_human), nsets = 14, order.by = "freq",
      text.scale = c(1.5, 1.5, 1.5, 1.2, 1.2, 1.5))
dev.off()
pdf("plots/upset_sequin_ONT.pdf", height = 5, width = 12)
upset(fromList(DE_sequin[1:7]), nsets = 7, order.by = "freq",
      text.scale = c(1.5, 1.5, 1.5, 1.2, 1.2, 1.5))
dev.off()
pdf("plots/upset_sequin_Illumina.pdf", height = 5, width = 12)
upset(fromList(DE_sequin[8:14]), nsets = 7, order.by = "freq",
      text.scale = c(1.5, 1.5, 1.5, 1.2, 1.2, 1.5))
dev.off()
pdf("plots/upset_sequin.pdf", height = 7, width = 15)
upset(fromList(DE_sequin), nsets = 14, order.by = "freq",
      text.scale = c(1.5, 1.5, 1.5, 1.2, 1.2, 1.5))
dev.off()
# calculate and plot sequins FDR and TPR
anno <- read.table("/wehisan/home/allstaff/d/dong.x/annotation/sequins/rnasequin_isoforms_2.4.tsv", header = TRUE, stringsAsFactors = FALSE)
anno$logFC <- log(anno$MIX_B / anno$MIX_A)
trueDE <- anno$NAME[anno$logFC != 0]
FDR_sequin <- sapply(DE_sequin, function(x){
  FD = sum(!(x %in% trueDE))
  FD / length(x)
})
res_sequin$FDR <- FDR_sequin
pdf("plots/FDR_sequin.pdf", height = 5, width = 8)
ggplot(res_sequin, aes(x = libsize, y = FDR, colour = dataset, group = dataset)) +
  geom_line() +
  geom_point() +
  scale_colour_manual(values = c("#FCB344", "#438DAC")) +
  labs(x = "Number of reads", y = "False discovery rate", colour = "Dataset") +
  theme_bw() +
  theme(text = element_text(size = 20))
dev.off()
# sequins TPR
TPR_sequin <- sapply(DE_sequin, function(x){
  TP = sum(x %in% trueDE)
  TP / length(trueDE)
})
res_sequin$TPR <- TPR_sequin
pdf("plots/TPR_sequin.pdf", height = 5, width = 8)
ggplot(res_sequin, aes(x = libsize, y = TPR, colour = dataset, group = dataset)) +
  geom_line() +
  geom_point() +
  scale_colour_manual(values = c("#FCB344", "#438DAC")) +
  labs(x = "Number of reads", y = "True positive rate", colour = "Dataset") +
  theme_bw() +
  theme(text = element_text(size = 20))
dev.off()
# calculate recovery method 1: treat the union of DE from full dataset as ground truth
DE_full_all <- union(DE_human$Illumina_full, DE_human$ONT_full)
recovery_human <- sapply(DE_human, function(x){
  sum(x %in% DE_full_all) / length(DE_full_all)
})
res_human$recovery <- recovery_human
pdf("plots/recovery_union_human.pdf", height = 5, width = 8)
ggplot(res_human, aes(x = libsize, y = recovery, colour = dataset, group = dataset)) +
  geom_line() +
  geom_point() +
  scale_colour_manual(values = c("#FCB344", "#438DAC")) +
  labs(x = "Number of reads", y = "Recovery rate", colour = "Dataset") +
  theme_bw() +
  theme(text = element_text(size = 20))
dev.off()
# calculate recovery rate method 2: use the intersect of DE from full datasets as ground truth
DE_full_intersect <- intersect(DE_human$Illumina_full, DE_human$ONT_full)
recovery_human_2 <- sapply(DE_human, function(x){
  sum(x %in% DE_full_intersect) / length(DE_full_intersect)
})
res_human$recovery2 <- recovery_human_2
pdf("plots/recovery_intersect_human.pdf", height = 5, width = 8)
ggplot(res_human, aes(x = libsize, y = recovery2, colour = dataset, group = dataset)) +
  geom_line() +
  geom_point() +
  scale_colour_manual(values = c("#FCB344", "#438DAC")) +
  labs(x = "Number of reads", y = "Recovery rate", colour = "Dataset") +
  theme_bw() +
  theme(text = element_text(size = 20))
dev.off()
# calculate recovery rate method 3: for each dataset, compare to the DE from its full dataset
res_human$recovery3[res_human$dataset == "ONT"] <- sapply(DE_human[res_human$dataset == "ONT"], function(x){
  sum(x %in% DE_human$ONT_full) / length(DE_human$ONT_full)
})
res_human$recovery3[res_human$dataset == "Illumina"] <- sapply(DE_human[res_human$dataset == "Illumina"], function(x){
  sum(x %in% DE_human$Illumina_full) / length(DE_human$Illumina_full)
})
pdf("plots/recovery_by_dataset_human.pdf", height = 5, width = 8)
ggplot(res_human, aes(x = libsize, y = recovery3, colour = dataset, group = dataset)) +
  geom_line() +
  geom_point() +
  scale_colour_manual(values = c("#FCB344", "#438DAC")) +
  labs(x = "Number of reads", y = "Recovery rate", colour = "Dataset") +
  theme_bw() +
  theme(text = element_text(size = 20))
dev.off()
# calculate inconsistency and make bar plot of consistency and inconsistency
res_human$consistent[res_human$dataset == "ONT"] <- sapply(DE_human[res_human$dataset == "ONT"], function(x){
  sum(x %in% DE_human$ONT_full)
})
res_human$consistent[res_human$dataset == "Illumina"] <- sapply(DE_human[res_human$dataset == "Illumina"], function(x){
  sum(x %in% DE_human$Illumina_full)
})
res_human$inconsistent[res_human$dataset == "ONT"] <- sapply(DE_human[res_human$dataset == "ONT"], function(x){
  sum(!(x %in% DE_human$ONT_full))
})
res_human$inconsistent[res_human$dataset == "Illumina"] <- sapply(DE_human[res_human$dataset == "Illumina"], function(x){
  sum(!(x %in% DE_human$Illumina_full))
})
res_human_DE <- data.table::melt(res_human,
                                 id = c(2, 3), measure = c(7, 8))
pdf("plots/DE_bar.pdf", height = 5, width = 10)
ggplot(res_human_DE, aes(x = libsize, y = ifelse(variable=="consistent", value, -value), 
                         fill = variable)) +
  geom_bar(stat = "identity") +
  facet_grid(cols = vars(dataset)) +
  scale_y_continuous(labels = abs, expand = expansion(mult = c(0.1, 0.1))) +
  labs(x = "Number of reads", y = "Number of DE transcripts",
       fill = "Category") +
  geom_text(aes(label = value,
            vjust = ifelse(variable=="consistent", -0.5, 1))) +
  scale_fill_manual(values = c("powderblue", "rosybrown2")) +
  theme_bw() +
  theme(text = element_text(size = 20))
dev.off()