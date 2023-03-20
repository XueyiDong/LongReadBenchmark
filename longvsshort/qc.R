library(ggplot2)
library(edgeR)
library(MetBrewer)
library(scales)
library(RColorBrewer)
library(tidyverse)
library(cowplot)
library(viridis)

DIR="/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark"
# load DGE lists
# s <- catchSalmon(file.path(DIR, "ONT/salmon_bs", list.files(file.path(DIR, "/ONT/salmon_bs"))))
# dge <- DGEList(counts=s$counts, genes=s$annotation)
# s.short <- catchSalmon(file.path(DIR, "illumina/salmon_bs", list.files(file.path(DIR, "illumina/salmon_bs"))))
# dge.short <- DGEList(counts = s.short$counts, genes = s.short$annotation)

dge <- readRDS("dge.rds")
dge.short <- readRDS("dge.short.rds")

# read num plot----
# organize read num stat 
read.stat <- data.frame(
  sample = rep(c(paste("H1975", 1:3, sep = "-"),
              paste("HCC827", c(1, 2, 5), sep = "-")), 2),
  raw_reads = c(46809934, 40748569, 42074278, 51044877, 48219467, 39145316,
               31073747, 34059318, 27974052, 31329213, 134000622, 30331521),
  mapped_reads = c(41054341, 35939259, 36437436, 46878103, 40629774, 36962945,
                   28999906, 31495290, 26126733, 29046678, 125097661, 26408681),
  # read_counts = c(dge$samples$lib.size[1:6], dge.short$samples$lib.size[1:6]),
  dataset = rep(c("ONT", "Illumina"), c(6,6))
)

read.stat <- data.table::melt(read.stat, id.vars = c("sample", "dataset"))
# read num plot
pdf("plots/readNum.pdf", height = 6, width = 6)
ggplot(read.stat, aes(x=variable, y=value, fill=sample, label = value))+
  geom_bar(stat="identity") +
  facet_grid(cols=vars(dataset)) +
  geom_text(aes(label = label_number_si()(value)), size = 3, colour = "white", position = position_stack(vjust = 0.5)) +
  ylab("Number") +
  xlab("Category") +
  theme_bw() +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6))+
  scale_fill_manual(values = met.brewer("Troy", 6))
dev.off()

# biotype ONT----
# gene biotype info
dge.human <- dge[grep("^ENST", rownames(dge)), ]
txid <- strsplit2(rownames(dge.human$counts), "\\|")[,1]
txid <- strsplit2(txid, "\\.")[,1]

library(AnnotationHub)
ah <- AnnotationHub()
EnsDb.Hsapiens.v104 <- query(ah, c("EnsDb", "Homo Sapiens", 104))[[1]]
biotype<- mapIds(
  x = EnsDb.Hsapiens.v104,
  # NOTE: Need to remove gene version number prior to lookup.
  keys = txid,
  keytype = "TXID",
  column = "TXBIOTYPE")
dge.human$genes$biotype <- biotype
# deal with biotype
# http://asia.ensembl.org/info/genome/genebuild/biotypes.html
dge.human$genes$biotype[grepl("pseudogene$", dge.human$genes$biotype)] <- "pseudogene"
dge.human$genes$biotype[grepl("^TR", dge.human$genes$biotype)] <- "IG_or_TR_gene"
dge.human$genes$biotype[grepl("^IG", dge.human$genes$biotype)] <- "IG_or_TR_gene"
dge.human$genes$biotype[dge.human$genes$biotype %in% c("miRNA", "misc_RNA", 
                                                       "piRNA", "rRNA", "siRNA",
                                                       "snRNA", "snoRNA", "scaRNA",
                                                       "tRNA", "vault_RNA", "scRNA",
                                                       "sRNA", "Mt_rRNA", "Mt_tRNA",
                                                       "ribozyme"
                                                       )] <- "ncRNA"
saveRDS(dge.human$genes, "txInfo.long.RDS")

# for each sample
biotype_sum <- sapply(1:6, function(x){
  typesum = aggregate(dge.human$counts[,x], by=list(dge.human$genes$biotype), FUN=sum, simplify=TRUE)
  return(typesum)
}, simplify=FALSE)
biotype_sum <- do.call("rbind", biotype_sum)
biotype_sum$sample <- rep(c("H1975-1", "H1975-2", "H1975-3", "HCC827-1", "HCC827-2", "HCC827-5"), rep(10, 6))
colnames(biotype_sum) <- c("biotype", "total_count", "sample")
# order the bars
ord = aggregate(biotype_sum$total_count, by = list(biotype_sum$biotype), FUN = sum, simplify = TRUE)
ord = ord[order(ord$x), ]

pdf("plots/biotype.pdf", height = 5, width = 8)
ggplot(biotype_sum, aes(x=sample, y=total_count, fill=factor(biotype, levels=ord$Group.1))) +
  geom_bar(stat="identity", position = "fill") +
  theme_bw() +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  labs(fill = "Transcript biotype", x = "Sample", y = "Proportion of count")
dev.off()

# calculate proportion for ONT----
biotype_sum$proportion <- sapply(1:nrow(biotype_sum), function(x){
  sample.sum = sum(biotype_sum$total_count[biotype_sum$sample == biotype_sum$sample[x]])
  return(biotype_sum$total_count[x] / sample.sum)
}, simplify = TRUE)

#biotype Illumina---- 
dge.short.human <- dge[grep("^ENST", rownames(dge.short)), ]
txid <- strsplit2(rownames(dge.short.human$counts), "\\|")[,1]
txid <- strsplit2(txid, "\\.")[,1]

biotype<- mapIds(
  x = EnsDb.Hsapiens.v104,
  # NOTE: Need to remove gene version number prior to lookup.
  keys = txid,
  keytype = "TXID",
  column = "TXBIOTYPE")
dge.short.human$genes$biotype <- biotype
saveRDS(biotype, "biotype.RDS")
# deal with biotype
dge.short.human$genes$biotype[grepl("pseudogene$", dge.short.human$genes$biotype)] <- "pseudogene"
dge.short.human$genes$biotype[grepl("^TR", dge.short.human$genes$biotype)] <- "IG_or_TR_gene"
dge.short.human$genes$biotype[grepl("^IG", dge.short.human$genes$biotype)] <- "IG_or_TR_gene"
dge.short.human$genes$biotype[dge.short.human$genes$biotype %in% c("miRNA", "misc_RNA", 
                                                       "piRNA", "rRNA", "siRNA",
                                                       "snRNA", "snoRNA", "scaRNA",
                                                       "tRNA", "vault_RNA", "scRNA",
                                                       "sRNA", "Mt_rRNA", "Mt_tRNA",
                                                       "ribozyme"
)] <- "ncRNA"
saveRDS(dge.short.human$genes, "txInfo.short.RDS")

# for each sample
biotype_sum.short <- sapply(1:6, function(x){
  typesum = aggregate(dge.short.human$counts[,x], by=list(dge.short.human$genes$biotype), FUN=sum, simplify=TRUE)
  return(typesum)
}, simplify=FALSE)
biotype_sum.short <- do.call("rbind", biotype_sum.short)

biotype_sum.short$sample <- rep(c("H1975-1", "H1975-2", "H1975-3", "HCC827-1", "HCC827-2", "HCC827-5"), rep(10, 6))
# biotype_sum.short$sample <- rep(paste0("barcode0", 1:6), rep(10, 6))
colnames(biotype_sum.short) <- c("biotype", "total_count", "sample")
# order the bars
ord = aggregate(biotype_sum$total_count, by = list(biotype_sum$biotype), FUN = sum, simplify = TRUE)
ord = ord[order(ord$x), ]

pdf("plots/biotype_short.pdf", height = 5, width = 8)
ggplot(biotype_sum.short, aes(x=sample, y=total_count, fill=factor(biotype, levels=ord$Group.1))) +
  geom_bar(stat="identity", position = "fill") +
  theme_bw() +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  labs(fill = "Transcript biotype", x = "Sample", y = "Proportion of count")
dev.off()

# calculate proportion for Illumina----
biotype_sum.short$proportion <- sapply(1:nrow(biotype_sum.short), function(x){
  sample.sum = sum(biotype_sum.short$total_count[biotype_sum.short$sample == biotype_sum.short$sample[x]])
  return(biotype_sum.short$total_count[x] / sample.sum)
}, simplify = TRUE)

#biotype long and short---- 
biotype_sum.all <- rbind(biotype_sum, biotype_sum.short)
biotype_sum.all$dataset <- rep(c("ONT", "Illumina"), c(nrow(biotype_sum), nrow(biotype_sum.short)))

pdf("plots/biotype_all.pdf", height = 4, width = 8)
ggplot(biotype_sum.all, aes(x=sample, y=total_count, fill=factor(biotype, levels=ord$Group.1))) +
  geom_bar(stat="identity", position = "fill") +
  facet_grid(cols=vars(dataset)) +
  theme_bw() +
  theme(text = element_text(size = 20), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        ) +
  scale_fill_brewer(palette = "Set3") +
  labs(fill = "Transcript biotype", x = "Sample", y = "Proportion of count")
dev.off()

# count distribution by biotype-----
dge.human$genes$counts <- rowSums(dge.human$counts[, 1:6])
dge.short.human$genes$counts <- rowSums(dge.short.human$counts[, 1:6])
dge.human$genes$dataset <- "ONT"
dge.short.human$genes$dataset <- "Illumina"
genes <- rbind(dge.human$genes, dge.short.human$genes)
genes <- genes[genes$counts>0, ]
genes <- na.omit(genes)
# library(ggridges)
pdf("plots/biotype_count.pdf", height = 4, width = 8)
ggplot(genes, aes(x=biotype, y=counts, fill=factor(biotype, levels=ord$Group.1))) +
  geom_violin() +
  geom_boxplot(width = 0.2, outlier.colour = NA, alpha = 0) +
  facet_grid(cols = vars(dataset)) +
  labs(x = "Transcript biotype", y = "Total count") +
  scale_fill_brewer(palette = "Set3") +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
# long vs short quantification------
# long CPM vs short TPM
# filter
dge.pure <- dge[,1:6]
dge.pure$samples$group <- rep(c("H1975", "HCC827"), c(3, 3))
dge.pure <- dge.pure[filterByExpr(dge.pure),] %>% calcNormFactors
cpm.long <- cpm(dge.pure)

dge.short.pure <- dge.short[,1:6]
dge.short.pure$samples$group <- rep(c("H1975", "HCC827"), c(3, 3))
dge.short.pure <- dge.short.pure[filterByExpr(dge.short.pure), ]%>% calcNormFactors
tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}
tpm.short <- tpm3(dge.short.pure$counts, dge.short.pure$genes$Length)
m <- match(rownames(dge.short.pure), rownames(dge.pure))
quant <- data.frame(
  TPM_short = c(log(rowMeans(tpm.short[, 1:3]) + 0.5),
                log(rowMeans(tpm.short[, 4:6]) + 0.5)),
  CPM_long = c(log(rowMeans(cpm.long[, 1:3]) + 0.5)[m],
               log(rowMeans(cpm.long[, 4:6]) + 0.5)[m]),
  group = rep(c("H1975", "HCC827"), rep(nrow(tpm.short), 2))
)
quant <- na.omit(quant)
# correlation
cor(quant$TPM_short, quant$CPM_long, method = "spearman")
cor(quant$TPM_short[quant$group=="H1975"], quant$CPM_long[quant$group=="H1975"], method = "spearman")
cor(quant$TPM_short[quant$group=="HCC827"], quant$CPM_long[quant$group=="HCC827"], method = "spearman")
pdf("plots/longVsShortQuant.pdf", height = 8, width = 9)
ggplot(quant, aes(x = CPM_long, y = TPM_short))+
  stat_binhex(bins=100) +
  scale_fill_viridis(trans = "log10", option = "A")+
  annotate(geom="text", x=3, y=12,
           label=paste0("Spearman's rho =", round(cor(quant$TPM_short, quant$CPM_long, method = "spearman"), 3)),
           size=7)+
  labs(x = expression("log"[2]*"CPM ONT read counts"),
       y = expression("log"[2]*"TPM Illumina read counts")
  ) +
  theme_bw() +
  theme(text=element_text(size = 20))
dev.off()

# quantification correlation matrix heatmap by biotype-------
## use raw count -----
m <- match(rownames(dge.short), rownames(dge))
table(is.na(m))
dge.all <- DGEList(counts = cbind(dge$counts[m, 1:6], dge.short$counts[, 1:6]))
dge.all$samples$group <- rep(c("H1975_long", "HCC827_long", "H1975_short", "HCC827_short"), c(3, 3, 3, 3))
# check filtering strategies
filt <- filterByExpr(dge.all)
table(filt)
filt.long <- filterByExpr(dge[, 1:6], group = rep(c("H1975", "HCC827"), each = 3))
filt.short <- filterByExpr(dge.short[, 1:6], group = rep(c("H1975", "HCC827"), each = 3))
filt2 <- filt.long[m] | filt.short
table(filt2)
table(filt == filt2)
# use the union of transcripts passed filtering in long and short read data
dge.all <- dge.all[filt2,]
dge.all$genes <- dge.short$genes[match(rownames(dge.all), rownames(dge.short)),]
cormat <- cor(dge.all$counts, method = "spearman")
library(pheatmap)
anno <- data.frame(
  cell_type = rep(rep(c("H1975", "HCC827"), c(3, 3)), 2),
  dataset = rep(c("ONT", "Illumina"), c(6, 6))
)
anno_colours = list(
  cell_type = c(H1975 = "#8b3a2b", HCC827 = "#235070"),
  dataset = c(ONT = "#438DAC", Illumina = "#FCB344")
)
rownames(anno) <- rownames(cormat)
pheatmap(cormat,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = anno,
         annotation_row = anno,
         annotation_colors = anno_colours,
         scale = "none",
         display_numbers = TRUE
         )
## use long cpm vs short tpm----
### all----
cpm <- cpm(dge.all[, 1:6])
tpm <- tpm3(dge.all$counts[, 7:12], dge.short$genes$Length[match(rownames(dge.all), rownames(dge.short))])
quant.all <- cbind(cpm, tpm)
cormat2 <- cor(quant.all, method = "spearman")
breaksList <- seq(0, 1, by = 0.01)
pdf("plots/corHeatmap.pdf", height = 8, width = 9)
pheatmap(cormat2,
         color = colorRampPalette(brewer.pal(n = 7, name = "PuBuGn"))(length(breaksList)),
         breaks = breaksList,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = anno,
         annotation_row = anno,
         annotation_colors = anno_colours,
         scale = "none",
         display_numbers = TRUE,
         number_color = "black",
         fontsize = 16,
         fontsize_number = 12
)
dev.off()
### coding----
# filter by biotype
tx.sel <- names(biotype)[biotype %in% c("protein_coding", "lncRNA")]
dge.pure.sel <- dge.pure[substr(rownames(dge.pure), 1, 15) %in% tx.sel, ]
dge.short.pure.sel <- dge.short.pure[substr(rownames(dge.short.pure), 1, 15) %in% tx.sel, ]
cpm.sel <- cpm(dge.pure.sel)
tpm.sel <- tpm3(dge.short.pure.sel$counts, dge.short.pure.sel$genes$Length)
quant.all.sel <- cbind(cpm.sel[match(rownames(tpm.sel), rownames(cpm.sel)), ], tpm.sel) %>% na.omit
cormat3 <- cor(quant.all.sel, method = "spearman")
pdf("plots/corHeatmapCoding.pdf", height = 7, width = 8)
pheatmap(cormat3,
         color = colorRampPalette(brewer.pal(n = 7, name = "PuBuGn"))(length(breaksList)),
         breaks = breaksList,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = anno,
         annotation_row = anno,
         annotation_colors = anno_colours,
         scale = "none",
         display_numbers = TRUE,
         number_color = "black",
         fontsize = 16,
         fontsize_number = 12
)
dev.off()

### other----
tx.sel <- names(biotype)[!(biotype %in% c("protein_coding", "lncRNA"))]
dge.pure.sel <- dge.pure[substr(rownames(dge.pure), 1, 15) %in% tx.sel, ]
dge.short.pure.sel <- dge.short.pure[substr(rownames(dge.short.pure), 1, 15) %in% tx.sel, ]
cpm.sel <- cpm(dge.pure.sel)
tpm.sel <- tpm3(dge.short.pure.sel$counts, dge.short.pure.sel$genes$Length)
quant.all.sel <- cbind(cpm.sel[match(rownames(tpm.sel), rownames(cpm.sel)), ], tpm.sel) %>% na.omit
cormat3 <- cor(quant.all.sel, method = "spearman")
pdf("plots/corHeatmapOther.pdf", height = 8, width = 9)
pheatmap(cormat3,
         color = colorRampPalette(brewer.pal(n = 7, name = "PuBuGn"))(length(breaksList)),
         breaks = breaksList,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = anno,
         annotation_row = anno,
         annotation_colors = anno_colours,
         scale = "none",
         display_numbers = TRUE,
         number_color = "black",
         fontsize = 16,
         fontsize_number = 12
)
dev.off()

### sequins----
cormat4 <- cor(quant.all[grepl("^R", rownames(quant.all)),], method = "spearman")
pdf("plots/corHeatmapSequin.pdf", height = 8, width = 9)
pheatmap(cormat4,
         color = colorRampPalette(brewer.pal(n = 7, name = "PuBuGn"))(length(breaksList)),
         breaks = breaksList,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = anno,
         annotation_row = anno,
         annotation_colors = anno_colours,
         scale = "none",
         display_numbers = TRUE,
         number_color = "black",
         fontsize = 16,
         fontsize_number = 12
)
dev.off()

### by expression quantile -----
# High: top 1/3. Median: 1/3~2/3. Low: bottom 1/3.
# Use total count rather than total TPM to calculate cutoffs.
#### low ----
cutoffs = quantile(rowSums(dge.all$counts), probs = c(1/3, 2/3))
summary(rowSums(dge.all$counts))
cutoffs
cormat.low <- cor(quant.all[rowSums(dge.all$counts) <= cutoffs[1], ], method = "spearman")
pdf("plots/corHeatmapLowExpr.pdf", height = 8, width = 9)
pheatmap(cormat.low,
         color = colorRampPalette(brewer.pal(n = 7, name = "PuBuGn"))(length(breaksList)),
         breaks = breaksList,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = anno,
         annotation_row = anno,
         annotation_colors = anno_colours,
         scale = "none",
         display_numbers = TRUE,
         number_color = "black",
         fontsize = 16,
         fontsize_number = 12
)
dev.off()
#### Median ----
cormat.med <- cor(quant.all[rowSums(dge.all$counts) > cutoffs[1] & rowSums(quant.all) <= cutoffs[2], ], method = "spearman")
pdf("plots/corHeatmapMedExpr.pdf", height = 8, width = 9)
pheatmap(cormat.med,
         color = colorRampPalette(brewer.pal(n = 7, name = "PuBuGn"))(length(breaksList)),
         breaks = breaksList,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = anno,
         annotation_row = anno,
         annotation_colors = anno_colours,
         scale = "none",
         display_numbers = TRUE,
         number_color = "black",
         fontsize = 16,
         fontsize_number = 12
)
dev.off()
#### High ----
cormat.high <- cor(quant.all[rowSums(dge.all$counts) > cutoffs[2], ], method = "spearman")
pdf("plots/corHeatmapHighExpr.pdf", height = 8, width = 9)
pheatmap(cormat.high,
         color = colorRampPalette(brewer.pal(n = 7, name = "PuBuGn"))(length(breaksList)),
         breaks = breaksList,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = anno,
         annotation_row = anno,
         annotation_colors = anno_colours,
         scale = "none",
         display_numbers = TRUE,
         number_color = "black",
         fontsize = 16,
         fontsize_number = 12
)
dev.off()
#### the formation of biotype in each expression category ----
biotype.expr <- data.frame(
 transcript = names(biotype),
 biotype = biotype
)
biotype.expr$biotype[grepl("pseudogene$", biotype.expr$biotype)] <- "pseudogene"
biotype.expr$biotype[grepl("^TR", biotype.expr$biotype)] <- "IG_or_TR_gene"
biotype.expr$biotype[grepl("^IG", biotype.expr$biotype)] <- "IG_or_TR_gene"
biotype.expr$biotype[biotype.expr$biotype %in% c("miRNA", "misc_RNA", 
                                                                   "piRNA", "rRNA", "siRNA",
                                                                   "snRNA", "snoRNA", "scaRNA",
                                                                   "tRNA", "vault_RNA", "scRNA",
                                                                   "sRNA", "Mt_rRNA", "Mt_tRNA",
                                                                   "ribozyme"
)] <- "ncRNA"
biotype.expr$total_count <- rowSums(dge.all$counts)[match(biotype.expr$transcript, substr(rownames(dge.all), 1, 15))]
biotype.expr <- na.omit(biotype.expr)
biotype.expr$expression_level[biotype.expr$total_count <= cutoffs[1]] <- "low"
biotype.expr$expression_level[biotype.expr$total_count > cutoffs[1] & biotype.expr$total_count <= cutoffs[2]] <- "medium"
biotype.expr$expression_level[biotype.expr$total_count > cutoffs[2]] <- "high"
biotype.expr$expression_level <- factor(biotype.expr$expression_level, levels = c("low", "medium", "high"))
pdf("plots/biotypeByExprGrp.pdf", height = 4, width = 8)
ggplot(biotype.expr, aes(fill = factor(biotype, levels=ord$Group.1), x = expression_level)) +
  geom_bar(position = "stack") +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "Expression level group", y = "Number of transcripts", fill = "Transcript biotype") +
  theme_bw() +
  theme(text = element_text(size = 16))
dev.off()
### by number of tx-------
# ### 1----
# cormat.tx1 <- cor(quant.all[dge.all$genes$nTranscript==1,], method = "spearman")
# pdf("plots/corHeatmaptx1.pdf", height = 8, width = 9)
# pheatmap(cormat.tx1,
#          color = colorRampPalette(brewer.pal(n = 7, name = "PuBuGn"))(length(breaksList)),
#          breaks = breaksList,
#          cluster_cols = FALSE,
#          cluster_rows = FALSE,
#          show_colnames = FALSE,
#          show_rownames = FALSE,
#          annotation_col = anno,
#          annotation_row = anno,
#          annotation_colors = anno_colours,
#          scale = "none",
#          display_numbers = TRUE,
#          number_color = "black",
#          fontsize = 16,
#          fontsize_number = 12
# )
# dev.off()
# 
# ### 2 - 5----
# cormat.tx2_5 <- cor(quant.all[dge.all$genes$nTranscript>=2 & dge.all$genes$nTranscript <= 5,], method = "spearman")
# pdf("plots/corHeatmaptx2_5.pdf", height = 8, width = 9)
# pheatmap(cormat.tx2_5,
#          color = colorRampPalette(brewer.pal(n = 7, name = "PuBuGn"))(length(breaksList)),
#          breaks = breaksList,
#          cluster_rows = FALSE,
#          show_colnames = FALSE,
#          show_rownames = FALSE,
#          annotation_col = anno,
#          annotation_row = anno,
#          annotation_colors = anno_colours,
#          scale = "none",
#          display_numbers = TRUE,
#          number_color = "black",
#          fontsize = 16,
#          fontsize_number = 12
# )
# dev.off()
# 
# ### > 5----
# cormat.tx5 <- cor(quant.all[dge.all$genes$nTranscript > 5,], method = "spearman")
# pdf("plots/corHeatmaptx5.pdf", height = 8, width = 9)
# pheatmap(cormat.tx5,
#          color = colorRampPalette(brewer.pal(n = 7, name = "PuBuGn"))(length(breaksList)),
#          breaks = breaksList,
#          cluster_cols = FALSE,
#          cluster_rows = FALSE,
#          show_colnames = FALSE,
#          show_rownames = FALSE,
#          annotation_col = anno,
#          annotation_row = anno,
#          annotation_colors = anno_colours,
#          scale = "none",
#          display_numbers = TRUE,
#          number_color = "black",
#          fontsize = 16,
#          fontsize_number = 12
# )
# dev.off()

# length bias plot--------------------
dge.pure$genes$totalCount <- rowSums(dge.pure$counts)
lm.long <- lm(dge.pure$genes$totalCount~dge.pure$genes$Length)
summary(lm.long)
# Multiple R-squared:  0.0007356,	Adjusted R-squared:  0.0007191 
pdf("plots/longLenBias.pdf", height = 5, width = 8)
ggplot(dge.pure$genes, aes(x = Length, y = totalCount))+
  stat_binhex() +
  geom_smooth(formula = y~x, method="lm") +
  scale_fill_viridis(trans = "log10", option = "A")+
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  annotate(geom="text", x=max(dge.pure$genes$Length) * 0.1, y=max(dge.pure$genes$totalCount) * 0.8,
           label=paste0("Pearson's r=", round(cor(dge.pure$genes$Length, dge.pure$genes$totalCount), 3)),
           size=7)+
  labs(x = "Transcript length",
       y = "Total read count"
  ) +
  theme_bw() +
  theme(text=element_text(size = 20))
dev.off()

p <- lapply(1:6, function(x){
  dat <- data.frame(
    Length = dge.pure$genes$Length,
    Count = dge.pure$counts[,x]
  )
  p <- ggplot(dat, aes(x = Length, y = Count))+
    stat_binhex() +
    geom_smooth(formula = y~x, method="lm") +
    scale_fill_viridis(trans = "log10", option = "A")+
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    annotate(geom="text", x=max(dat$Length) * 0.05, y=max(dat$Count) * 0.8,
             label=paste0("Pearson's r=", round(cor(dat$Length, dat$Count), 3)),
             size=7)+
    labs(x = "Transcript length",
         y = "Total read count"
    ) +
    theme_bw() +
    theme(text=element_text(size = 20), axis.text.x = element_text(angle = 30, hjust = 1))
  return(p)
})
pdf("plots/longLenBiasSamp.pdf", height = 9, width = 16)
plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]],
          labels = paste(rep(c("H1975", "HCC827"), c(3, 3)), c(1, 2, 3, 1, 2, 5), sep = "-"))
dev.off()

dge.short.pure$genes$totalCount <- rowSums(dge.short.pure$counts)
lm.short <- lm(dge.short.pure$genes$totalCount~dge.short.pure$genes$Length)
summary(lm.short)
# Multiple R-squared:  3.275e-06,	Adjusted R-squared:  -2.818e-05
pdf("plots/shortLenBias.pdf", height = 5, width = 8)
ggplot(dge.short.pure$genes, aes(x = Length, y = totalCount))+
  stat_binhex() +
  geom_smooth(formula = y~x, method="lm") +
  scale_fill_viridis(trans = "log10", option = "A")+
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  annotate(geom="text", x=max(dge.short.pure$genes$Length) * 0.1, y=max(dge.short.pure$genes$totalCount) * 0.8,
           label=paste0("Pearson's r=", round(cor(dge.short.pure$genes$Length, dge.short.pure$genes$totalCount), 3)),
           size=7)+
  labs(x = "Transcript length",
       y = "Total read count"
  ) +
  theme_bw() +
  theme(text=element_text(size = 20))
dev.off()

p <- lapply(1:6, function(x){
  dat <- data.frame(
    Length = dge.short.pure$genes$Length,
    Count = dge.short.pure$counts[,x]
  )
  p <- ggplot(dat, aes(x = Length, y = Count))+
    stat_binhex() +
    geom_smooth(formula = y~x, method="lm") +
    scale_fill_viridis(trans = "log10", option = "A")+
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    annotate(geom="text", x=max(dat$Length) * 0.05, y=max(dat$Count) * 0.8,
             label=paste0("Pearson's r=", round(cor(dat$Length, dat$Count), 3)),
             size=7)+
    labs(x = "Transcript length",
         y = "Total read count"
    ) +
    theme_bw() +
    theme(text=element_text(size = 20), axis.text.x = element_text(angle = 30, hjust = 1))
  return(p)
})
pdf("plots/ShortLenBiasSamp.pdf", height = 9, width = 16)
plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]],
          labels = paste(rep(c("H1975", "HCC827"), c(3, 3)), c(1, 2, 3, 1, 2, 5), sep = "-"))
dev.off()

# BCV plots ---------
## Calculate human and sequin together -----
### ONT -------
dge.pure <- dge[, 1:6]
dge.pure$samples$group = rep(c("H1975", "HCC827"), each = 3)
dge.pure <- dge.pure[filterByExpr(dge.pure),]
dge.pure <- dge.pure %>% 
  calcNormFactors() %>% 
  estimateDisp()

plotBCV(dge.pure)

dispersion <- data.frame(
  AveLogCPM = dge.pure$AveLogCPM,
  common.dispersion = dge.pure$common.dispersion,
  trended.dispersion = dge.pure$trended.dispersion,
  tagwise.dispersion = dge.pure$tagwise.dispersion
)
rownames(dispersion) <- rownames(dge.pure)
dispersion <- dispersion[order(dispersion$AveLogCPM), ]

plot.bcv.long <- ggplot() +
  stat_binhex(data = dispersion[grep("^ENST", rownames(dispersion)), ],
              aes(x = AveLogCPM, y = sqrt(tagwise.dispersion)), bins = 100) +
  scale_fill_viridis(option = "magma") +
  geom_line(data = dispersion[grep("^ENST", rownames(dispersion)), ],
            aes(x = AveLogCPM, y = sqrt(common.dispersion), colour = "Common")) +
  geom_line(data = dispersion[grep("^ENST", rownames(dispersion)), ],
            aes(x = AveLogCPM, y = sqrt(trended.dispersion), colour = "Trend")) +
  geom_point(data = dispersion[grep("^R", rownames(dispersion)), ],
             aes(x = AveLogCPM, y = sqrt(tagwise.dispersion),
                 colour = "Tagwise_sequins")) +
  scale_colour_manual(name = "BCV",
                      breaks = c("Common", "Trend", "Tagwise_sequins"),
                      values = c("Common" = "red", "Trend" = "blue", "Tagwise_sequins" = "green")) +
  scale_y_continuous(limits = c(0, 0.75)) +
  theme_bw() +
  labs(x = "Average log CPM", y = "Bioligical coefficient of variation",
       fill = "Tagwise BCV\ncount",
       title = "ONT long-read BCV") +
  theme(text = element_text(size = 16))

### Illumina ---------
dge.short.pure <- dge.short[, 1:6]
dge.short.pure$samples$group = rep(c("H1975", "HCC827"), each = 3)
dge.short.pure <- dge.short.pure[filterByExpr(dge.short.pure),]
dge.short.pure <- dge.short.pure %>% 
  calcNormFactors() %>% 
  estimateDisp()

plotBCV(dge.short.pure)

dispersion.short <- data.frame(
  AveLogCPM = dge.short.pure$AveLogCPM,
  common.dispersion = dge.short.pure$common.dispersion,
  trended.dispersion = dge.short.pure$trended.dispersion,
  tagwise.dispersion = dge.short.pure$tagwise.dispersion
)
rownames(dispersion.short) <- rownames(dge.short.pure)
dispersion.short <- dispersion.short[order(dispersion.short$AveLogCPM), ]

plot.bcv.short <- ggplot() +
  stat_binhex(data = dispersion.short[grep("^ENST", rownames(dispersion.short)), ],
              aes(x = AveLogCPM, y = sqrt(tagwise.dispersion)), bins = 100) +
  scale_fill_viridis(option = "magma") +
  geom_line(data = dispersion.short[grep("^ENST", rownames(dispersion.short)), ],
            aes(x = AveLogCPM, y = sqrt(common.dispersion), colour = "Common")) +
  geom_line(data = dispersion.short[grep("^ENST", rownames(dispersion.short)), ],
            aes(x = AveLogCPM, y = sqrt(trended.dispersion), colour = "Trend")) +
  geom_point(data = dispersion.short[grep("^R", rownames(dispersion.short)), ],
             aes(x = AveLogCPM, y = sqrt(tagwise.dispersion),
                 colour = "Tagwise_sequins")) +
  scale_colour_manual(name = "BCV",
                      breaks = c("Common", "Trend", "Tagwise_sequins"),
                      values = c("Common" = "red", "Trend" = "blue", "Tagwise_sequins" = "green")) +
  scale_y_continuous(limits = c(0, 0.75)) +
  theme_bw() +
  labs(x = "Average log CPM", y = "Bioligical coefficient of variation",
       fill = "Tagwise BCV\ncount",
       title = "Illumina short-read BCV") +
  theme(text = element_text(size = 16))

### PDF ----
pdf("plots/bcv.pdf", height = 5, width = 15)
cowplot::plot_grid(plot.bcv.long, plot.bcv.short)
dev.off()

## Calculate human and sequins separately ----
### ONT ----
dge.pure.human <- dge.pure[grepl("^ENST", rownames(dge.pure)), ] %>% 
  calcNormFactors() %>% 
  estimateDisp()
dispersion.human <- data.frame(
  AveLogCPM = dge.pure.human$AveLogCPM,
  common.dispersion = dge.pure.human$common.dispersion,
  trended.dispersion = dge.pure.human$trended.dispersion,
  tagwise.dispersion = dge.pure.human$tagwise.dispersion
)
rownames(dispersion.human) <- rownames(dge.pure.human)
dispersion.human <- dispersion.human[order(dispersion.human$AveLogCPM), ]

dge.pure.sequin <- dge.pure[grepl("^R", rownames(dge.pure)), ] %>% 
  calcNormFactors() %>% 
  estimateDisp()
dispersion.sequin <- data.frame(
  AveLogCPM = dge.pure.sequin$AveLogCPM,
  common.dispersion = dge.pure.sequin$common.dispersion,
  trended.dispersion = dge.pure.sequin$trended.dispersion,
  tagwise.dispersion = dge.pure.sequin$tagwise.dispersion
)
rownames(dispersion.sequin) <- rownames(dge.pure.sequin)
dispersion.sequin <- dispersion.sequin[order(dispersion.sequin$AveLogCPM), ]

plot.bcv.long.sep <- 
  ggplot() +
  stat_binhex(data =dispersion.human,
              aes(x = AveLogCPM, y = sqrt(tagwise.dispersion)), bins = 100) +
  scale_fill_viridis(option = "magma") +
  geom_line(data = dispersion.human,
            aes(x = AveLogCPM, y = sqrt(common.dispersion), 
                colour = "Common", linetype = "Human")) +
  geom_line(data = dispersion.human,
            aes(x = AveLogCPM, y = sqrt(trended.dispersion),
                colour = "Trend", linetype = "Human")) +
  geom_point(data = dispersion.sequin,
             aes(x = AveLogCPM, y = sqrt(tagwise.dispersion),
                 colour = "Tagwise_sequins")) +
  geom_line(data = dispersion.sequin,
            aes(x = AveLogCPM, y = sqrt(common.dispersion), 
                colour = "Common", linetype = "Sequin")) +
  geom_line(data = dispersion.sequin,
            aes(x = AveLogCPM, y = sqrt(trended.dispersion),
                colour = "Trend", linetype = "Sequin")) +
  scale_colour_manual(name = "BCV",
                      breaks = c("Common", "Trend", "Tagwise_sequins"),
                      values = c("Common" = "red", "Trend" = "blue", "Tagwise_sequins" = "seagreen")) +
  scale_linetype_manual(name = "Source",
                        breaks = c("Human", "Sequin"),
                        values = c("Human" = 1, "Sequin" = 2)) +
  scale_y_continuous(limits = range(c(sqrt(dispersion.human$tagwise.dispersion), sqrt(dispersion.sequin$tagwise.dispersion)))) +
  theme_bw() +
  labs(x = "Average log CPM", y = "Bioligical coefficient of variation",
       fill = "Tagwise BCV\ncount",
       title = "ONT long-read BCV") +
  theme(text = element_text(size = 16))

### Illumina -----
dge.short.pure.human <- dge.short.pure[grepl("^ENST", rownames(dge.short.pure)), ] %>% 
  calcNormFactors() %>% 
  estimateDisp()
dispersion.short.human <- data.frame(
  AveLogCPM = dge.short.pure.human$AveLogCPM,
  common.dispersion = dge.short.pure.human$common.dispersion,
  trended.dispersion = dge.short.pure.human$trended.dispersion,
  tagwise.dispersion = dge.short.pure.human$tagwise.dispersion
)
rownames(dispersion.short.human) <- rownames(dge.short.pure.human)
dispersion.short.human <- dispersion.short.human[order(dispersion.short.human$AveLogCPM), ]

dge.short.pure.sequin <- dge.short.pure[grepl("^R", rownames(dge.short.pure)), ] %>% 
  calcNormFactors() %>% 
  estimateDisp()
dispersion.short.sequin <- data.frame(
  AveLogCPM = dge.short.pure.sequin$AveLogCPM,
  common.dispersion = dge.short.pure.sequin$common.dispersion,
  trended.dispersion = dge.short.pure.sequin$trended.dispersion,
  tagwise.dispersion = dge.short.pure.sequin$tagwise.dispersion
)
rownames(dispersion.short.sequin) <- rownames(dge.short.pure.sequin)
dispersion.short.sequin <- dispersion.short.sequin[order(dispersion.short.sequin$AveLogCPM), ]

plot.bcv.short.sep <- 
  ggplot() +
  stat_binhex(data =dispersion.short.human,
              aes(x = AveLogCPM, y = sqrt(tagwise.dispersion)), bins = 100) +
  scale_fill_viridis(option = "magma") +
  geom_line(data = dispersion.short.human,
            aes(x = AveLogCPM, y = sqrt(common.dispersion), 
                colour = "Common", linetype = "Human")) +
  geom_line(data = dispersion.short.human,
            aes(x = AveLogCPM, y = sqrt(trended.dispersion), 
                colour = "Trend", linetype = "Human")) +
  geom_point(data = dispersion.short.sequin,
             aes(x = AveLogCPM, y = sqrt(tagwise.dispersion),
                 colour = "Tagwise_sequins")) +
  geom_line(data = dispersion.short.sequin,
            aes(x = AveLogCPM, y = sqrt(common.dispersion), 
                colour = "Common", linetype = "Sequin")) +
  geom_line(data = dispersion.short.sequin,
            aes(x = AveLogCPM, y = sqrt(trended.dispersion), 
                colour = "Trend", linetype = "Sequin")) +
  scale_colour_manual(name = "BCV",
                      breaks = c("Common", "Trend", "Tagwise_sequins"),
                      values = c("Common" = "red", "Trend" = "blue", "Tagwise_sequins" = "seagreen")) +
  scale_y_continuous(limits = range(c(sqrt(dispersion.short.human$tagwise.dispersion), sqrt(dispersion.short.sequin$tagwise.dispersion)))) +
  scale_linetype_manual(name = "Source",
                        breaks = c("Human", "Sequin"),
                        values = c("Human" = 1, "Sequin" = 2)) +
  theme_bw() +
  labs(x = "Average log CPM", y = "Bioligical coefficient of variation",
       fill = "Tagwise BCV\ncount",
       title = "Illumina short-read BCV") +
  theme(text = element_text(size = 16))

## Combined plot ----
# ggplot
pdf("plots/bcv2.pdf", height = 5, width = 16)
cowplot::plot_grid(plot.bcv.long.sep, 
                   plot.bcv.short.sep,
                   labels = c("a", "b"),
                   label_size = 24
                   )
dev.off()
# base R plot
pdf("plots/bcv_base.pdf", height = 5, width = 16)
par(mfrow = c(1, 2))
plotBCV(dge.pure.human,
        ylim = range(c(sqrt(dispersion.human$tagwise.dispersion), sqrt(dispersion.sequin$tagwise.dispersion))),
        main = "ONT long-read BCV")
points(dispersion.sequin$AveLogCPM, sqrt(dispersion.sequin$common.dispersion),
       pch = 16, cex = .5, col = "seagreen")
abline(h = sqrt(dispersion.sequin$common.dispersion), col = "red", lwd = 2, lty = 2)
lines(dispersion.sequin$AveLogCPM, sqrt(dispersion.sequin$trended.dispersion), col = "blue", lwd = 2, lty = 2)
legend("topright", legend = c("Tagwise human", "Common human", "Trend human", 
                              "Tagwise sequin", "Common sequin", "Trend sequin"),
       lty = c(-1, 1, 1, -1, 2, 2),
       pch = c(16, -1, -1, 16, -1, -1),
       pt.cex = c(.2, 1, 1, .5, 1, 1),
       lwd = 2,
       col = c("black", "red", "blue", "seagreen", "red", "blue"),
       bg = "white")

plotBCV(dge.short.pure.human,
        ylim = range(c(sqrt(dispersion.short.human$tagwise.dispersion), sqrt(dispersion.short.sequin$tagwise.dispersion))),
        main = "Illumina short-read BCV")
points(dispersion.short.sequin$AveLogCPM, sqrt(dispersion.short.sequin$common.dispersion),
       pch = 16, cex = .5, col = "seagreen")
abline(h = sqrt(dispersion.short.sequin$common.dispersion), col = "red", lwd = 2, lty = 2)
lines(dispersion.short.sequin$AveLogCPM, sqrt(dispersion.short.sequin$trended.dispersion), col = "blue", lwd = 2, lty = 2)
legend("topright", legend = c("Tagwise human", "Common human", "Trend human", 
                              "Tagwise sequin", "Common sequin", "Trend sequin"),
       lty = c(-1, 1, 1, -1, 2, 2),
       pch = c(16, -1, -1, 16, -1, -1),
       pt.cex = c(.2, 1, 1, .5, 1, 1),
       lwd = 2,
       col = c("black", "red", "blue", "seagreen", "red", "blue"),
       bg = "white")
dev.off()
