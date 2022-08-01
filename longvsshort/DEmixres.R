library(ggplot2)
library(viridis)

DIR="/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark"

DE.human.ONT <- readRDS(file.path(DIR, "ONT/DTEmix/DE.human.RDS"))
DE.human.illumina <- readRDS(file.path(DIR, "illumina/DTEmix/DE.human.RDS"))
DE.sequin.ONT <- readRDS(file.path(DIR, "ONT/DTEmix/DE.sequin.RDS"))
DE.sequin.illumina <- readRDS(file.path(DIR, "illumina/DTEmix/DE.sequin.RDS"))
# names(DE.human.ONT) <- c("limma", "edgeR", "DESeq2", "EBSeq", "NOISeq")
# names(DE.human.illumina) <- c("limma", "edgeR", "DESeq2", "EBSeq", "NOISeq")
# names(DE.sequin.ONT) <- c("limma", "edgeR", "DESeq2", "EBSeq", "NOISeq")
# names(DE.sequin.illumina) <- c("limma", "edgeR", "DESeq2", "EBSeq", "NOISeq")

DEgenecomp <- data.frame(
  method = character(0),
  comparison = character(0),
  category = character(0),
  number = numeric(0)
)
for(i in 1:5){
  for(x in 1:4){
    ill = DE.human.illumina[[i]][[x]]
    ont = DE.human.ONT[[i]][[x]]
    int = intersect(ill, ont)
    common = length(int)
    Illumina_only = length(ill) - common
    ONT_only = length(ont) - common
    tmp = data.frame(
      method = rep(c("limma", "edgeR", "DESeq2", "EBSeq", "NOISeq")[i], 3),
      comparison = rep(c("100vs000", "075vs025", "050vs025", "075vs050")[x], 3),
      category = c("both", "Illumina_only", "ONT_only"),
      number = c(common, Illumina_only, ONT_only)
    )
    DEgenecomp <- rbind(DEgenecomp, tmp)
  }
}
DEgenecomp$comparison <- factor(DEgenecomp$comparison, levels=c("100vs000", "075vs025", "050vs025", "075vs050"))

# bar plot
pdf("plots/DTE/DTEbarHuman.pdf", height = 5, width = 18)
ggplot(DEgenecomp, aes(x=method, y=number, fill=comparison))+
  geom_bar(stat="identity", position="dodge") +
  facet_grid(cols=vars(category)) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  scale_fill_manual(values=c("#CB664F", "#D3B1A7", "#8DBFA3", "#DEBA47"))
dev.off()

# venn diagrams for long vs short of each method 100 vs 000
library(VennDiagram)
for(i in 1:5){
  venn.diagram(
    x = list(DE.human.illumina[[i]][[1]], DE.human.ONT[[i]][[1]]),
    category.names = c("Illumina", "ONT"),
    filename = paste("plots/DTEvenn_", c("limma", "edgeR", "DESeq2", "EBSeq", "NOISeq")[i],".png", sep=""),
    output = FALSE,
    imagetype="png" ,
    height = 800 , 
    width = 800, 
    resolution = 300,
    compression = "lzw",
    cat.dist = c(0.05,0.05),
    cat.default.pos = "outer",
    fontfamily = "sans",
    cat.pos = c(-150, 150),
    margin = 0.06,
    ext.text=TRUE,
    na="remove",
    main = c("limma", "edgeR", "DESeq2", "EBSeq", "NOISeq")[i],
    col = c("#FCB344", "#438DAC"),
    fill = c("#FCB344", "#438DAC"),
    cat.col = c("#FCB344", "#438DAC"),
    alpha = 0.4
  )
}

# upset plot for long vs short 100 vs 000
# human
DE.human.illumina.100vs000 <- lapply(DE.human.illumina, function(x){
  return(x[[1]])
})
names(DE.human.illumina.100vs000) <- paste(c("limma", "edgeR", "DESeq2", "EBSeq", "NOISeq"), "Illumina", sep = "_")
DE.human.ONT.100vs000 <- lapply(DE.human.ONT, function(x){
  return(x[[1]])
})
names(DE.human.ONT.100vs000) <- paste(c("limma", "edgeR", "DESeq2", "EBSeq", "NOISeq"), "ONT", sep = "_")
saveRDS(DE.human.illumina.100vs000, "DE.human.illumina.100vs000.RDS")
saveRDS(DE.human.ONT.100vs000, "DE.human.ONT.100vs000.RDS")

library(UpSetR)
pdf("plots/DTE/DTEhumanUpset.pdf", height = 5, width = 11)
upset(fromList(append(DE.human.illumina.100vs000, DE.human.ONT.100vs000)), 
      nsets=10, nintersects = 25, order.by = "freq",
      text.scale = c(1.5, 1.5, 1.5, 1.2, 1.2, 1.5),
      sets.bar.color = rep(c("#D96A70", "#476937",  "#D5A2CB", "#708FA6", "#9FC675"), 2)[order(sapply(append(DE.human.illumina.100vs000, DE.human.ONT.100vs000), length, simplify = T), decreasing = TRUE)])
      # sets.bar.color = rep(c("#D96A70", "#476937",  "#D5A2CB", "#708FA6", "#9FC675"), 2))
dev.off()
pdf("plots/DTE/DTEhumanUpsetONT.pdf", height = 5, width = 8)
upset(fromList(DE.human.ONT.100vs000), order.by = "freq",
      sets.bar.color = c("#D96A70", "#476937",  "#D5A2CB", "#708FA6", "#9FC675")[order(sapply(DE.human.ONT.100vs000, length, simplify = T), decreasing = T)])
dev.off()
pdf("plots/DTE/DTEhumanUpsetIllumina.pdf", height = 5, width = 8)
upset(fromList(DE.human.illumina.100vs000), order.by = "freq",
      sets.bar.color = c("#D96A70", "#476937",  "#D5A2CB", "#708FA6", "#9FC675")[order(sapply(DE.human.illumina.100vs000, length, simplify = T), decreasing = T)])
dev.off()




# sequin
DE.sequin.illumina.100vs000 <- lapply(DE.sequin.illumina, function(x){
  return(x[[1]])
})
names(DE.sequin.illumina.100vs000) <- paste(c("limma", "edgeR", "DESeq2", "EBSeq", "NOISeq"), "Illumina", sep = "_")
DE.sequin.ONT.100vs000 <- lapply(DE.sequin.ONT, function(x){
  return(x[[1]])
})
names(DE.sequin.ONT.100vs000) <- paste(c("limma", "edgeR", "DESeq2", "EBSeq", "NOISeq"), "ONT", sep = "_")

saveRDS(DE.sequin.illumina.100vs000, "DE.sequin.illumina.100vs000.RDS")
saveRDS(DE.sequin.ONT.100vs000, "DE.sequin.ONT.100vs000.RDS")

pdf("plots/DTE/DTEsequinUpset.pdf", height = 5, width = 8)
upset(fromList(append(DE.sequin.illumina.100vs000, DE.sequin.ONT.100vs000)), 
      nsets=10, nintersects = 25, order.by = "freq",
      text.scale = c(1.5, 1.5, 1.5, 1.2, 1.2, 1.5),
      sets.bar.color = rep(c("#D96A70", "#476937",  "#D5A2CB", "#708FA6", "#9FC675"), 2)[order(sapply(append(DE.sequin.illumina.100vs000, DE.sequin.ONT.100vs000), length, simplify = T), decreasing = TRUE)])
dev.off()

# human and sequin
DE.illumina.100vs000 <- DE.human.illumina.100vs000
DE.ONT.100vs000 <- DE.human.ONT.100vs000
for(i in 1:5){
  DE.illumina.100vs000[[i]] <- c(DE.illumina.100vs000[[i]], DE.sequin.illumina.100vs000[[i]])
  DE.ONT.100vs000[[i]] <- c(DE.ONT.100vs000[[i]], DE.sequin.ONT.100vs000[[i]])
}
pdf("plots/DTE/DTEUpset.pdf", height = 5, width = 11)
upset(fromList(append(DE.illumina.100vs000, DE.ONT.100vs000)), 
      nsets=10, nintersects = 25, order.by = "freq",
      text.scale = c(1.5, 1.5, 1.5, 1.2, 1.2, 1.5),
      sets.bar.color = rep(c("#D96A70", "#476937",  "#D5A2CB", "#708FA6", "#9FC675"), 2)[order(sapply(append(DE.illumina.100vs000, DE.ONT.100vs000), length, simplify = T), decreasing = TRUE)])
dev.off()

# tx tested by all methods only
tx.human.ONT <- readRDS(file.path(DIR, "ONT/DTEmix/tx.DE.human.RDS"))
tx.sequin.ONT <- readRDS(file.path(DIR, "ONT/DTEmix/tx.DE.sequin.RDS"))
tx.human.illumina <- readRDS(file.path(DIR, "illumina/DTEmix/tx.DE.human.RDS"))
tx.sequin.illumina <- readRDS(file.path(DIR, "illumina/DTEmix/tx.DE.sequin.RDS"))
tx.human.ONT.100vs000 <- lapply(tx.human.ONT, function(x){
  return(x[[1]])
})
tx.sequin.ONT.100vs000 <- lapply(tx.sequin.ONT, function(x){
  return(x[[1]])
})
tx.human.illumina.100vs000 <- lapply(tx.human.illumina, function(x){
  return(x[[1]])
})
tx.sequin.illumina.100vs000 <- lapply(tx.sequin.illumina, function(x){
  return(x[[1]])
})
tx.human <- Reduce(intersect, append(tx.human.ONT.100vs000, tx.human.illumina.100vs000))
tx.sequin <- Reduce(intersect, append(tx.sequin.ONT.100vs000, tx.sequin.illumina.100vs000))
DE.illumina.100vs000.filt <- lapply(DE.illumina.100vs000, function(x){
  return(x[x %in% c(tx.human, tx.sequin)])
})
DE.ONT.100vs000.filt <- lapply(DE.ONT.100vs000, function(x){
  return(x[x %in% c(tx.human, tx.sequin)])
})
pdf("plots/DTE/DTEUpsetFilt.pdf", height = 7, width = 16)
upset(fromList(append(DE.illumina.100vs000.filt, DE.ONT.100vs000.filt)), 
      nsets=10, nintersects = 25, order.by = "freq",
      text.scale = c(2, 2, 2, 1.75, 1.8, 2),
      sets.bar.color = rep(c("#D96A70", "#476937",  "#D5A2CB", "#708FA6", "#9FC675"), 2)[order(sapply(append(DE.illumina.100vs000.filt, DE.ONT.100vs000.filt), length, simplify = T), decreasing = TRUE)])
dev.off()

# long vs short t
tt.human.ONT <- read.table(file.path(DIR, "ONT/DTEmix/topTableHumanc100vs0.tsv"), sep = "\t", header = T)
tt.human.illumina <- read.table(file.path(DIR, "illumina/DTEmix/topTableHumanc100vs0.tsv"), sep = "\t", header = T)
tt.sequin.ONT <- read.table(file.path(DIR, "ONT/DTEmix/topTableSequinc100vs0.tsv"), sep = "\t", header = T)
tt.sequin.illumina <- read.table(file.path(DIR, "illumina/DTEmix/topTableSequinc100vs0.tsv"), sep = "\t", header = T)
m <- match(tt.human.illumina$TXNAME, tt.human.ONT$TXNAME)
m2 <- match(tt.sequin.illumina$TXNAME, tt.sequin.ONT$TXNAME)

t <- data.frame(
  TXNAME = c(tt.human.ONT$TXNAME[m], tt.sequin.ONT$TXNAME[m2]),
  t.long = c(tt.human.ONT$t[m], tt.sequin.ONT$t[m2]),
  t.short = c(tt.human.illumina$t, tt.sequin.illumina$t),
  source = rep(c("human", "sequin"), c(nrow(tt.human.illumina), nrow(tt.sequin.illumina)))
)
t$z.long <- limma::zscoreT(t$t.long, df=4)
t$z.short <- limma::zscoreT(t$t.short, df = 4)
t <- na.omit(t)
# add biotype information
t$biotype <- biotype[match(substr(t$TXNAME, 1, 15), names(biotype))]
t.filt <- t[t$biotype %in% c("protein_coding", "lncRNA"),]
t.filt2 <- t[!(t$biotype %in% c("protein_coding", "lncRNA")),]
pdf("plots/DTE/t.pdf", height = 5, width = 8)
ggplot(t, aes(x=t.long, y=t.short)) +
  stat_binhex(bins=100) +
  geom_smooth(method='lm', formula= y~x) +
  theme_bw() +
  labs(x="ONT t-statistic", y = "Illumina t-statistic", fill = "Density:\nnumber of \ntranscripts") +
  annotate(geom="text", x = max(t$t.long) * 0.6, y=max(t$t.short) * 0.95, 
           label=paste0("Pearson's r=", round(cor(t$t.long, t$t.short), 3)), size=7) +
  scale_fill_viridis(direction = 1, option="A", trans = "log10") +
  theme(text=element_text(size = 20))
ggplot(t.filt, aes(x=t.long, y=t.short)) +
  stat_binhex(bins=100) +
  geom_smooth(method='lm', formula= y~x) +
  theme_bw() +
  labs(x="ONT t-statistic", y = "Illumina t-statistic", fill = "Density:\nnumber of \ntranscripts",
       title = "Protein coding and lncRNA") +
  annotate(geom="text", x = max(t.filt$t.long) * 0.6, y=max(t.filt$t.short) * 0.95, 
           label=paste0("Pearson's r=", round(cor(t.filt$t.long, t.filt$t.short), 3)), size=7) +
  scale_fill_viridis(direction = 1, option="A", trans = "log10") +
  theme(text=element_text(size = 20))
ggplot(t.filt2, aes(x=t.long, y=t.short)) +
  stat_binhex(bins=100) +
  geom_smooth(method='lm', formula= y~x) +
  theme_bw() +
  labs(x="ONT t-statistic", y = "Illumina t-statistic", fill = "Density:\nnumber of \ntranscripts",
       title = "Other biotypes") +
  annotate(geom="text", x = max(t.filt2$t.long) * 0.6, y=max(t.filt2$t.short) * 0.95, 
           label=paste0("Pearson's r=", round(cor(t.filt2$t.long, t.filt2$t.short), 3)), size=7) +
  scale_fill_viridis(direction = 1, option="A", trans = "log10") +
  theme(text=element_text(size = 20))
dev.off()

pdf("plots/DTE/z.pdf", height = 5, width = 8)
ggplot(t, aes(x=z.long, y=z.short)) +
  stat_binhex(bins=100) +
  geom_smooth(method='lm', formula= y~x) +
  theme_bw() +
  labs(x="ONT z-score of t-statistic", y = "Illumina z-score of t-statistic", fill = "Density:\nnumber of \ntranscripts") +
  annotate(geom="text", x = max(t$z.long) * 0.6, y = max(t$z.short) * 0.95, 
           label=paste0("Pearson's r=", round(cor(t$z.long, t$z.short), 3)), size=7) +
  scale_fill_viridis(direction = 1, option="A", trans = "log10") +
  theme(text=element_text(size = 20))
ggplot(t.filt, aes(x=z.long, y=z.short)) +
  stat_binhex(bins=100) +
  geom_smooth(method='lm', formula= y~x) +
  theme_bw() +
  labs(x="ONT z-score of t-statistic", y = "Illumina z-score of t-statistic", fill = "Density:\nnumber of \ntranscripts",
       title = "Protein coding and lncRNA") +
  annotate(geom="text", x = max(t.filt$z.long) * 0.6, y = max(t.filt$z.short) * 0.95, 
           label=paste0("Pearson's r=", round(cor(t.filt$z.long, t.filt$z.short), 3)), size=7) +
  scale_fill_viridis(direction = 1, option="A", trans = "log10") +
  theme(text=element_text(size = 20))
ggplot(t.filt2, aes(x=z.long, y=z.short)) +
  stat_binhex(bins=100) +
  geom_smooth(method='lm', formula= y~x) +
  theme_bw() +
  labs(x="ONT z-score of t-statistic", y = "Illumina z-score of t-statistic", fill = "Density:\nnumber of \ntranscripts",
       title = "Other biotypes") +
  annotate(geom="text", x = max(t.filt2$z.long) * 0.6, y = max(t.filt2$z.short) * 0.95, 
           label=paste0("Pearson's r=", round(cor(t.filt2$z.long, t.filt2$z.short), 3)), size=7) +
  scale_fill_viridis(direction = 1, option="A", trans = "log10") +
  theme(text=element_text(size = 20))
dev.off()

# biotype and length of DTEs
txInfo.long <- readRDS("txInfo.long.RDS")
txInfo.short <- readRDS("txInfo.short.RDS")
DE.human.100vs000 <- data.frame(
  tx = unlist(append(DE.human.ONT.100vs000, DE.human.illumina.100vs000), use.names = TRUE),
  method = rep(rep(c("limma", "edgeR", "DESeq2", "EBSeq", "NOISeq"), 2),
               sapply(append(DE.human.ONT.100vs000, DE.human.illumina.100vs000), length, simplify = TRUE)),
  dataset = rep(rep(c("ONT", "Illumina"), c(5, 5)),
                sapply(append(DE.human.ONT.100vs000, DE.human.illumina.100vs000), length, simplify = TRUE))
)
DE.human.100vs000$biotype <- c(txInfo.long$biotype[match(DE.human.100vs000$tx[DE.human.100vs000$dataset=="ONT"],
                                                         strsplit2(rownames(txInfo.long), "\\|")[,1])],
                               txInfo.short$biotype[match(DE.human.100vs000$tx[DE.human.100vs000$dataset=="Illumina"],
                                                          strsplit2(rownames(txInfo.short), "\\|")[,1])]
)
DE.human.100vs000$length <- c(txInfo.long$Length[match(DE.human.100vs000$tx[DE.human.100vs000$dataset=="ONT"],
                                                       strsplit2(rownames(txInfo.long), "\\|")[,1])],
                              txInfo.short$Length[match(DE.human.100vs000$tx[DE.human.100vs000$dataset=="Illumina"],
                                                        strsplit2(rownames(txInfo.short), "\\|")[,1])]
)
DE.human.100vs000 <- na.omit(DE.human.100vs000)
col <- RColorBrewer::brewer.pal(10, "Set3")
ord <- readRDS("ord.RDS")
pdf("plots/DTE/DTEbiotype.pdf", height = 5, width = 8)
ggplot(DE.human.100vs000, aes(x = method, fill=factor(biotype, levels=ord$Group.1)))+
  geom_bar(position = "fill")+
  facet_grid(cols=vars(dataset)) +
  theme_bw() +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette="Set3") +
  # scale_fill_manual(values = col[-1]) +
  labs(fill = "Transcript biotype", x = "Method", y = "Proportion of DTE transcripts")
dev.off()

pdf("plots/DTE/DTEbiotypeCount.pdf", height = 5, width = 8)
ggplot(DE.human.100vs000, aes(x = method, fill=factor(biotype, levels=ord$Group.1)))+
  geom_bar(position = "stack")+
  facet_grid(cols=vars(dataset)) +
  theme_bw() +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette="Set3") +
  # scale_fill_manual(values = col[-1]) +
  labs(fill = "Transcript biotype", x = "Method", y = "Number of DTE transcripts")
dev.off()

# explore biotype of long- or short-read only DTE
category <- readRDS("DTEcategory.RDS")
category2 <- readRDS("DTEcategory.short.RDS")
category <- rbind(category, category2)
DE.human.100vs000$category <- category$category[match(DE.human.100vs000$tx, rownames(category))]

pdf("plots/DTE/DTECategoryBiotype.pdf", height = 5, width = 8)
ggplot(na.omit(DE.human.100vs000[!duplicated(DE.human.100vs000$tx),]), aes(x = category, fill=factor(biotype, levels=ord$Group.1)))+
  geom_bar(position = "fill")+
  theme_bw() +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette="Set3") +
  # scale_fill_manual(values = col[-1]) +
  labs(fill = "Transcript biotype", x = "Category", y = "Number of DTE transcripts")
dev.off()

library(ggridges)
pdf("plots/DTE/DTElength.pdf", height = 5, width = 8)
ggplot(DE.human.100vs000, aes(x = length, y=method, fill=method)) +
  geom_density_ridges(alpha = .7) +
  scale_fill_manual(values = c( "#D5A2CB", "#708FA6", "#476937", "#D96A70", "#9FC675"))+
  scale_x_continuous(trans = "log10") +
  labs(x = "Annotated transcript length") +
  facet_grid(rows=vars(dataset)) +
  theme_bw() +
  theme(text = element_text(size = 20), legend.position = "NA")
dev.off()

