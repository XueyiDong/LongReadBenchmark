DTU.gene.human.ONT <- readRDS("../ONT/DTUmix/DTU.gene.human.RDS")
DTU.gene.sequin.ONT <- readRDS("../ONT/DTUmix/DTU.gene.sequin.RDS")
DTU.tx.human.ONT <- readRDS("../ONT/DTUmix/DTU.tx.human.RDS")
DTU.tx.sequin.ONT <- readRDS("../ONT/DTUmix/DTU.tx.sequin.RDS")
DTU.gene.human.illumina <- readRDS("../illumina/DTUmix/DTU.gene.human.RDS")
DTU.gene.sequin.illumina <- readRDS("../illumina/DTUmix/DTU.gene.sequin.RDS")
DTU.tx.human.illumina <- readRDS("../illumina/DTUmix/DTU.tx.human.RDS")
DTU.tx.sequin.illumina <- readRDS("../illumina/DTUmix/DTU.tx.sequin.RDS")

DTUgenecomp <- data.frame(
  method = character(0),
  comparison = character(0),
  category = character(0),
  number = numeric(0)
)
for(i in 1:5){
  for(x in 1:4){
    ill = DTU.gene.human.illumina[[i]][[x]]
    ont = DTU.gene.human.ONT[[i]][[x]]
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
    DTUgenecomp <- rbind(DTUgenecomp, tmp)
  }
}
DTUgenecomp$comparison <- factor(DTUgenecomp$comparison, levels=c("100vs000", "075vs025", "050vs025", "075vs050"))
# bar plot
pdf("plots/DTUbarGeneHuman.pdf", height = 5, width = 18)
ggplot(DTUgenecomp, aes(x=method, y=number, fill=comparison))+
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
    x = list(DTU.gene.human.illumina[[i]][[1]], DTU.gene.human.ONT[[i]][[1]]),
    category.names = c("Illumina", "ONT"),
    filename = paste("plots/DTUgeneVenn_", c("DEXSeq", "DRIMSeq", "edgeR", "limma", "satuRn")[i],".png", sep=""),
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
DTU.gene.human.illumina.100vs000 <- lapply(DTU.gene.human.illumina, function(x){
  return(x[[1]])
})
names(DTU.gene.human.illumina.100vs000) <- paste(c("DEXSeq", "DRIMSeq", "edgeR", "limma", "satuRn"), "Illumina", sep = "_")
DTU.gene.human.ONT.100vs000 <- lapply(DTU.gene.human.ONT, function(x){
  return(x[[1]])
})
names(DTU.gene.human.ONT.100vs000) <- paste(c("DEXSeq", "DRIMSeq", "edgeR", "limma", "satuRn"), "ONT", sep = "_")
library(UpSetR)
pdf("plots/DTUgenehumanUpset.pdf", height = 5, width = 11)
upset(fromList(append(DTU.gene.human.illumina.100vs000, DTU.gene.human.ONT.100vs000)), 
      nsets=10, nintersects = 25, order.by = "freq",
      text.scale = c(1.5, 1.5, 1.5, 1.2, 1.2, 1.5),
      sets.bar.color = rep(c("#ECD98B", "#AAAAC2",  "#03875C", "#9A4C43", "#4E3227"), 2)[order(sapply(append(DTU.gene.human.illumina.100vs000, DTU.gene.human.ONT.100vs000), length, simplify = T), decreasing = T)])
dev.off()

# sequin
DTU.gene.sequin.illumina.100vs000 <- lapply(DTU.gene.sequin.illumina, function(x){
  return(x[[1]])
})
names(DTU.gene.sequin.illumina.100vs000) <- paste(c("DEXSeq", "DRIMSeq", "edgeR", "limma", "satuRn"), "Illumina", sep = "_")
DTU.gene.sequin.ONT.100vs000 <- lapply(DTU.gene.sequin.ONT, function(x){
  return(x[[1]])
})
names(DTU.gene.sequin.ONT.100vs000) <- paste(c("DEXSeq", "DRIMSeq", "edgeR", "limma", "satuRn"), "ONT", sep = "_")
pdf("plots/DTUgenesequinUpset.pdf", height = 5, width = 11)
upset(fromList(append(DTU.gene.sequin.illumina.100vs000, DTU.gene.sequin.ONT.100vs000)), 
      nsets=10, nintersects = 25, order.by = "freq",
      text.scale = c(1.5, 1.5, 1.5, 1.2, 1.2, 1.5),
      sets.bar.color = c("#ECD98B", "#AAAAC2",  "#03875C", "#9A4C43", "#4E3227")[c(2, 1, 5, 4, 3, 2, 5, 1, 4, 3)])
dev.off()