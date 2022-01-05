library(ggplot2)
library(viridis)

DE.human.ONT <- readRDS("../ONT/DTEmix/DE.human.RDS")
DE.human.illumina <- readRDS("../illumina/DTEmix/DE.human.RDS")
DE.sequin.ONT <- readRDS("../ONT/DTEmix/DE.sequin.RDS")
DE.sequin.illumina <- readRDS("../illumina/DTEmix/DE.sequin.RDS")
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
pdf("plots/DTEbarHuman.pdf", height = 5, width = 18)
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
pdf("plots/DTEhumanUpset.pdf", height = 5, width = 11)
upset(fromList(append(DE.human.illumina.100vs000, DE.human.ONT.100vs000)), 
      nsets=10, nintersects = 25, order.by = "freq",
      text.scale = c(1.5, 1.5, 1.5, 1.2, 1.2, 1.5),
      sets.bar.color = rep(c("#D96A70", "#476937",  "#D5A2CB", "#708FA6", "#9FC675"), 2)[order(sapply(append(DE.human.illumina.100vs000, DE.human.ONT.100vs000), length, simplify = T), decreasing = TRUE)])
      # sets.bar.color = rep(c("#D96A70", "#476937",  "#D5A2CB", "#708FA6", "#9FC675"), 2))
dev.off()
pdf("plots/DTEhumanUpsetONT.pdf", height = 5, width = 8)
upset(fromList(DE.human.ONT.100vs000), order.by = "freq",
      sets.bar.color = c("#D96A70", "#476937",  "#D5A2CB", "#708FA6", "#9FC675")[order(sapply(DE.human.ONT.100vs000, length, simplify = T), decreasing = T)])
dev.off()
pdf("plots/DTEhumanUpsetIllumina.pdf", height = 5, width = 8)
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

pdf("plots/DTEsequinUpset.pdf", height = 5, width = 8)
upset(fromList(append(DE.sequin.illumina.100vs000, DE.sequin.ONT.100vs000)), nsets=10, nintersects = 25, order.by = "freq")
dev.off()

# long vs short t
tt.human.ONT <- read.table("../ONT/DTEmix/topTableHumanc100vs0.tsv", sep = "\t", header = T)
tt.human.illumina <- read.table("../illumina/DTEmix/topTableHumanc100vs0.tsv", sep = "\t", header = T)
tt.sequin.ONT <- read.table("../ONT/DTEmix/topTableSequinc100vs0.tsv", sep = "\t", header = T)
tt.sequin.illumina <- read.table("../illumina/DTEmix/topTableSequinc100vs0.tsv", sep = "\t", header = T)
m <- match(tt.human.illumina$TXNAME, tt.human.ONT$TXNAME)
m2 <- match(tt.sequin.illumina$TXNAME, tt.sequin.ONT$TXNAME)

t <- data.frame(
  t.long = c(tt.human.ONT$t[m], tt.sequin.ONT$t[m2]),
  t.short = c(tt.human.illumina$t, tt.sequin.illumina$t),
  source = rep(c("human", "sequin"), c(nrow(tt.human.illumina), nrow(tt.sequin.illumina)))
  
)
t$z.long <- limma::zscoreT(t$t.long, df=4)
t$z.short <- limma::zscoreT(t$t.short, df = 4)
t <- na.omit(t)

t.lm <- lm(t$t.short ~ t$t.long)
summary(t.lm)
z.lm <- lm(t$z.short ~ t$z.long)
summary(z.lm)

pdf("plots/t.pdf", height = 5, width = 8)
ggplot(t, aes(x=t.long, y=t.short)) +
  stat_binhex() +
  geom_smooth(method='lm', formula= y~x) +
  theme_bw() +
  labs(x="ONT t-statistic", y = "Illumina t-statistic", fill = "Density:\nnumber of \ntranscripts") +
  annotate(geom="text", x=40, y=100, label="Adj R2 = 0.56\np-value < 2.2e-16", size=6) +
  scale_fill_viridis(direction = -1, option="A", trans = "log10") +
  theme(text=element_text(size = 20)) 
dev.off()

pdf("plots/z.pdf", height = 5, width = 8)
ggplot(t, aes(x=z.long, y=z.short)) +
  stat_binhex() +
  geom_smooth(method='lm', formula= y~x) +
  theme_bw() +
  labs(x="ONT z-score", y = "Illumina z-score", fill = "Density:\nnumber of \ntranscripts") +
  annotate(geom="text", x = 3, y = 5, label="Adj R2 = 0.59\np-value < 2.2e-16", size=6) +
  scale_fill_viridis(direction = -1, option="A", trans = "log10") +
  theme(text=element_text(size = 20)) 
dev.off()