library(ggplot2)
library(UpSetR)

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


# Compare long and short human t-statistic
tt.human.long <- read.delim("../ONT/topTableHuman.tsv", sep= "\t", stringsAsFactors = FALSE)
tt.human.short <- read.delim("../illumina/topTableHuman.tsv", sep = "\t", stringsAsFactors = FALSE)
m <- match(tt.human.long$GeneID, tt.human.short$GeneID)
cor(tt.human.long$t, tt.human.short$t[m], use = "complete.obs")
pdf("plots/humant.pdf")
smoothScatter(tt.human.long$t, tt.human.short$t[m], 
              xlab = "long read t-statistics", ylab = "short read t-statistics")
abline(coef = c(0,1))
dev.off()
