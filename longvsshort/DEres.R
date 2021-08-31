library(ggplot2)
library(UpSetR)
library(viridis)
# library(wesanderson)

res.human.long <- readRDS("../ONT/DEres/res.human.RDS")
res.sequin.long <- readRDS("../ONT/DEres/res.sequin.RDS")
res.human.short <- readRDS("../illumina/DEres/res.human.RDS")
res.sequin.short <- readRDS("../illumina/DEres/res.sequin.RDS")

# Organize res.human
res.human <- data.frame()
for(i in 1:length(res.human.long)){
  res.human <- rbind(res.human, res.human.long[[i]])
}
for(i in 1:length(res.human.short)){
  res.human <- rbind(res.human, res.human.short[[i]])
}
res.human$method <- rep(c("DESeq2", "EBSeq", "edgeR", "limma", "NOISeq", "DESeq2", "EBSeq", "edgeR", "limma", "NOISeq"),
                        c(sapply(res.human.long, nrow, simplify=TRUE), sapply(res.human.short, nrow, simplify=TRUE)))
res.human$dataset <- rep(c("ONT", "Illumina"),
                         c(sum(sapply(res.human.long, nrow, simplify=TRUE)), sum(sapply(res.human.short, nrow, simplify=TRUE))))
res.human$method_dataset <- paste(res.human$method, res.human$dataset, sep="_")

# Human DE genes upset plot
DEgenes.human <- lapply(unique(res.human$method_dataset), function(x){
  res <- res.human[res.human$method_dataset==x, ]
  DEgenes <- res$Gene[res$FDR < 0.05]
  return(DEgenes)
})
names(DEgenes.human) <- unique(res.human$method_dataset)
pdf("plots/humanDeUpset.pdf", height = 5)
upset(fromList(DEgenes.human), nsets=10, nintersects = 25, order.by = "freq")
dev.off()

# Organize res.sequin
res.sequin <- data.frame()
for(i in 1:length(res.sequin.long)){
  res.sequin <- rbind(res.sequin, res.sequin.long[[i]])
}
for(i in 1:length(res.sequin.short)){
  res.sequin <- rbind(res.sequin, res.sequin.short[[i]])
}
res.sequin$method <- rep(c("DESeq2", "EBSeq", "edgeR", "limma", "NOISeq", "DESeq2", "EBSeq", "edgeR", "limma", "NOISeq"),
                         c(sapply(res.sequin.long, nrow, simplify=TRUE), sapply(res.sequin.short, nrow, simplify=TRUE)))
res.sequin$dataset <- rep(c("ONT", "Illumina"),
                          c(sum(sapply(res.sequin.long, nrow, simplify=TRUE)), sum(sapply(res.sequin.short, nrow, simplify=TRUE))))
res.sequin$method_dataset <- paste(res.sequin$method, res.sequin$dataset, sep="_")

# Sequin DE genes upset plot
DEgenes.sequin <- lapply(unique(res.sequin$method_dataset), function(x){
  res <- res.sequin[res.sequin$method_dataset==x, ]
  DEgenes <- res$Gene[res$FDR < 0.05]
  return(DEgenes)
})
names(DEgenes.sequin) <- unique(res.sequin$method_dataset)
pdf("plots/sequinDeUpset.pdf", height = 5)
upset(fromList(DEgenes.sequin), nsets=10, order.by = "freq")
dev.off()

# Sequin FDR and TPR
anno <- read.table("/wehisan/home/allstaff/d/dong.x/annotation/sequins/rnasequin_genes_2.4.tsv", header = TRUE, stringsAsFactors = FALSE)
anno$logFC <- log(anno$MIX_B / anno$MIX_A)
res.sequin$logFC_expected <- anno$logFC[match(res.sequin$Gene, anno$NAME)]
# FDR = false positive / predicted positive
FDR <- as.data.frame(t(sapply(unique(res.sequin$method_dataset), function(x){
  fdr = nrow(res.sequin[res.sequin$method_dataset == x & res.sequin$FDR < 0.05 & res.sequin$logFC_expected == 0,]) / nrow(res.sequin[res.sequin$method_dataset == x & res.sequin$FDR < 0.05,] )
  return(c(x, fdr))
}, simplify=TRUE)))
colnames(FDR) <- c("Method_dataset", "False discovery rate")
FDR$`False discovery rate` <- as.numeric(FDR$`False discovery rate`)
pdf("plots/DEsequinFDR.pdf", height = 4)
ggplot(FDR, aes(x=Method_dataset, y=`False discovery rate`, fill=Method_dataset)) +
  geom_bar(stat="identity") +
  theme_bw()+
  scale_fill_brewer(palette = "Paired") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
# TPR = true positive / positive
TPR <- as.data.frame(t(sapply(unique(res.sequin$method_dataset), function(x){
  tpr = nrow(res.sequin[res.sequin$method_dataset == x & res.sequin$FDR < 0.05 & res.sequin$logFC_expected != 0,]) / nrow(res.sequin[res.sequin$method_dataset == x & res.sequin$logFC_expected != 0,] )
  return(c(x, tpr))
}, simplify=TRUE)))
colnames(TPR) <- c("Method_dataset", "True positive rate")
TPR$`True positive rate` <- as.numeric(TPR$`True positive rate`)
pdf("plots/DEsequinTPR.pdf", height = 4)
ggplot(TPR, aes(x=Method_dataset, y=`True positive rate`, fill=Method_dataset)) +
  geom_bar(stat="identity") +
  theme_bw() +
  scale_fill_brewer(palette = "Paired") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# check false positive/false negative genes
FD <- res.sequin[res.sequin$FDR<0.05 & res.sequin$logFC_expected == 0, ]
FN <- res.sequin[res.sequin$FDR >= 0.05 & res.sequin$logFC_expected != 0, ]
# venn diagram 

# Plot FDR vs TPR
FDR$`True positive rate` <- TPR$`True positive rate`
FDR$Method <- strsplit2(FDR$Method_dataset, "_")[,1]
FDR$Dataset <- strsplit2(FDR$Method_dataset, "_")[,2]
FDR$Method <- factor(FDR$Method, levels = c("limma", "edgeR", "DESeq2", "EBSeq", "NOISeq"))
pdf("plots/DEsequinFDRvsTPR.pdf", height = 5, width = 8)
ggplot(FDR, aes(x=`False discovery rate`, y=`True positive rate`, colour=Method, shape=Dataset))+
  geom_jitter(size=7, alpha = 0.9) +
  scale_colour_manual(values = c("#D96A70", "#476937",  "#D5A2CB", "#708FA6", "#9FC675"))+
  theme_bw() +
  theme(text = element_text(size = 20))
dev.off()

# Compare long and short human t-statistic
tt.human.long <- read.delim("../ONT/topTableHuman.tsv", sep= "\t", stringsAsFactors = FALSE)
tt.human.short <- read.delim("../illumina/topTableHuman.tsv", sep = "\t", stringsAsFactors = FALSE)
m <- match(tt.human.long$GeneID, tt.human.short$GeneID)
cor(tt.human.long$t, tt.human.short$t[m], use = "complete.obs")
pdf("plots/humant.pdf")
smoothScatter(tt.human.long$t, tt.human.short$t[m], 
              xlab = "long read t-statistics", ylab = "short read t-statistics",
              cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, cex = 1.5)
# abline(coef = c(0,1))
dev.off()

tt.sequin.long <- read.delim("../ONT/topTableSequin.tsv", sep= "\t", stringsAsFactors = FALSE)
tt.sequin.short <- read.delim("../illumina/topTableSequin.tsv", sep= "\t", stringsAsFactors = FALSE)
m2 <- match(tt.sequin.long$GeneID, tt.sequin.short$GeneID)
t <- data.frame(
  t.long = c(tt.human.long$t, tt.sequin.long$t),
  t.short = c(tt.human.short$t[m], tt.sequin.short$t[m2]),
  source = rep(c("human", "sequin"), c(nrow(tt.human.long), nrow(tt.sequin.long)))
  
)

t.lm <- lm(t$t.short ~ t$t.long)
summary(t.lm)
pdf("plots/t.pdf", height = 5, width = 8)
ggplot(t, aes(x=t.long, y=t.short)) +
  stat_binhex() +
  geom_smooth(method='lm', formula= y~x) +
  theme_bw() +
  labs(x="ONT t-statistic", y = "Illumina t-statistic", fill = "Density:\nnumber of \ngenes") +
  annotate(geom="text", x=30, y=90, label="Adj R2 = 0.82\np-value < 2.2e-16", size=6) +
  scale_fill_viridis(direction = -1, option="A", trans = "log10") +
  theme(text=element_text(size = 20)) 
dev.off()
