DTU.gene.human.ONT <- readRDS("../ONT/DTUmix/DTU.gene.human.RDS")
DTU.gene.sequin.ONT <- readRDS("../ONT/DTUmix/DTU.gene.sequin.RDS")
DTU.tx.human.ONT <- readRDS("../ONT/DTUmix/DTU.tx.human.RDS")
DTU.tx.sequin.ONT <- readRDS("../ONT/DTUmix/DTU.tx.sequin.RDS")
DTU.gene.human.illumina <- readRDS("../illumina/DTUmix/DTU.gene.human.RDS")
DTU.gene.sequin.illumina <- readRDS("../illumina/DTUmix/DTU.gene.sequin.RDS")
DTU.tx.human.illumina <- readRDS("../illumina/DTUmix/DTU.tx.human.RDS")
DTU.tx.sequin.illumina <- readRDS("../illumina/DTUmix/DTU.tx.sequin.RDS")

DTUtxcomp <- data.frame(
  method = character(0),
  comparison = character(0),
  category = character(0),
  number = numeric(0)
)
for(i in 1:5){
  for(x in 1:4){
    ill = DTU.tx.human.illumina[[i]][[x]]
    ont = DTU.tx.human.ONT[[i]][[x]]
    int = intersect(ill, ont)
    common = length(int)
    Illumina_only = length(ill) - common
    ONT_only = length(ont) - common
    tmp = data.frame(
      method = rep(c("DEXSeq", "DRIMSeq", "edgeR", "limma", "satuRn")[i], 3),
      comparison = rep(c("100vs000", "075vs025", "050vs025", "075vs050")[x], 3),
      category = c("both", "Illumina_only", "ONT_only"),
      number = c(common, Illumina_only, ONT_only)
    )
    DTUtxcomp <- rbind(DTUtxcomp, tmp)
  }
}
DTUtxcomp$comparison <- factor(DTUtxcomp$comparison, levels=c("100vs000", "075vs025", "050vs025", "075vs050"))
# bar plot
pdf("plots/DTUbarGeneHuman.pdf", height = 5, width = 18)
ggplot(DTUtxcomp, aes(x=method, y=number, fill=comparison))+
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
    x = list(DTU.tx.human.illumina[[i]][[1]], DTU.tx.human.ONT[[i]][[1]]),
    category.names = c("Illumina", "ONT"),
    filename = paste("plots/DTUtxVenn_", c("DEXSeq", "DRIMSeq", "edgeR", "limma", "satuRn")[i],".png", sep=""),
    output = FALSE,
    imagetype="png" ,
    height = 800 , 
    width = 800, 
    resolution = 300,
    compression = "lzw",
    cat.dist = c(0.05,0.05),
    cat.default.pos = "outer",
    fontfamily = "sans",
    cat.pos = c(-30, 30),
    margin = 0.06,
    ext.text=TRUE,
    na="remove",
    main = c("DEXSeq", "DRIMSeq", "edgeR", "limma", "satuRn")[i],
    col = c("#FCB344", "#438DAC"),
    fill = c("#FCB344", "#438DAC"),
    cat.col = c("#FCB344", "#438DAC"),
    alpha = 0.4
  )
}

# upset plot for long vs short 100 vs 000
# human
DTU.tx.human.illumina.100vs000 <- lapply(DTU.tx.human.illumina, function(x){
  return(x[[1]])
})
names(DTU.tx.human.illumina.100vs000) <- paste(c("DEXSeq", "DRIMSeq", "edgeR", "limma", "satuRn"), "Illumina", sep = "_")
DTU.tx.human.ONT.100vs000 <- lapply(DTU.tx.human.ONT, function(x){
  return(x[[1]])
})
names(DTU.tx.human.ONT.100vs000) <- paste(c("DEXSeq", "DRIMSeq", "edgeR", "limma", "satuRn"), "ONT", sep = "_")
library(UpSetR)
pdf("plots/DTUtxhumanUpset.pdf", height = 5, width = 8)
upset(fromList(append(DTU.tx.human.illumina.100vs000, DTU.tx.human.ONT.100vs000)), 
      nsets=10, nintersects = 25, order.by = "freq",
      text.scale = c(1.5, 1.5, 1.5, 1.2, 1.2, 1.5),
      sets.bar.color = rep(c("#ECD98B", "#AAAAC2",  "#03875C", "#9A4C43", "#4E3227"), 2)[order(sapply(append(DTU.tx.human.illumina.100vs000, DTU.tx.human.ONT.100vs000), length, simplify = T), decreasing = T)])
dev.off()

# sequin
DTU.tx.sequin.illumina.100vs000 <- lapply(DTU.tx.sequin.illumina, function(x){
  return(x[[1]])
})
names(DTU.tx.sequin.illumina.100vs000) <- paste(c("DEXSeq", "DRIMSeq", "edgeR", "limma", "satuRn"), "Illumina", sep = "_")
DTU.tx.sequin.ONT.100vs000 <- lapply(DTU.tx.sequin.ONT, function(x){
  return(x[[1]])
})
names(DTU.tx.sequin.ONT.100vs000) <- paste(c("DEXSeq", "DRIMSeq", "edgeR", "limma", "satuRn"), "ONT", sep = "_")
pdf("plots/DTUtxsequinUpset.pdf", height = 5, width = 8)
upset(fromList(append(DTU.tx.sequin.illumina.100vs000, DTU.tx.sequin.ONT.100vs000)), 
      nsets=10, nintersects = 25, order.by = "freq",
      text.scale = c(1.5, 1.5, 1.5, 1.2, 1.2, 1.5),
      sets.bar.color = c("#ECD98B", "#AAAAC2",  "#03875C", "#9A4C43", "#4E3227")[c(2, 4, 1, 5, 5, 2, 3, 1, 4, 3)]
      )
dev.off()

# long vs short t
ts.human.ONT <- read.table("../ONT/DTUmix/topSpliceHumanTxc100vs0.tsv", sep = "\t", header = T)
ts.human.illumina <- read.table("../illumina/DTUmix/topSpliceHumanTxc100vs0.tsv", sep = "\t", header = T)
ts.sequin.ONT <- read.table("../ONT/DTUmix/topSpliceSequinTxc100vs0.tsv", sep = "\t", header = T)
ts.sequin.illumina <- read.table("../illumina/DTUmix/topSpliceSequinTxc100vs0.tsv", sep = "\t", header = T)
m <- match(ts.human.illumina$ExonID, ts.human.ONT$ExonID)
m2 <- match(ts.sequin.illumina$ExonID, ts.sequin.ONT$ExonID)

t <- data.frame(
  t.long = c(ts.human.ONT$t[m], ts.sequin.ONT$t[m2]),
  t.short = c(ts.human.illumina$t, ts.sequin.illumina$t),
  source = rep(c("human", "sequin"), c(nrow(ts.human.illumina), nrow(ts.sequin.illumina)))
  
)
t$z.long <- limma::zscoreT(t$t.long, df=4)
t$z.short <- limma::zscoreT(t$t.short, df = 4)
t <- na.omit(t)

t.lm <- lm(t$t.short ~ t$t.long)
summary(t.lm)
z.lm <- lm(t$z.short ~ t$z.long)
summary(z.lm)

pdf("plots/t_DTU.pdf", height = 5, width = 8)
ggplot(t, aes(x=t.long, y=t.short)) +
  stat_binhex() +
  geom_smooth(method='lm', formula= y~x) +
  theme_bw() +
  labs(x="ONT t-statistic", y = "Illumina t-statistic", fill = "Density:\nnumber of \ntranscripts") +
  annotate(geom="text", x=40, y=100, label="Adj R2 = 0.26\np-value < 2.2e-16", size=6) +
  scale_fill_viridis(direction = -1, option="A", trans = "log10") +
  theme(text=element_text(size = 20)) 
dev.off()

pdf("plots/z_DTU.pdf", height = 5, width = 8)
ggplot(t, aes(x=z.long, y=z.short)) +
  stat_binhex() +
  geom_smooth(method='lm', formula= y~x) +
  theme_bw() +
  labs(x="ONT z-score", y = "Illumina z-score", fill = "Density:\nnumber of \ntranscripts") +
  annotate(geom="text", x = 3, y = 5, label="Adj R2 = 0.23\np-value < 2.2e-16", size=6) +
  scale_fill_viridis(direction = -1, option="A", trans = "log10") +
  theme(text=element_text(size = 20)) 
dev.off()