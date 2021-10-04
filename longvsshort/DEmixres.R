library(ggplot2)

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
  theme(text = element_text(size = 20))
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
