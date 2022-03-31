library(ggplot2)
library(UpSetR)

DIR="/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark"

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
pdf("plots/DTU/DTUbarGeneHuman.pdf", height = 5, width = 18)
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
pdf("plots/DTU/DTUgenehumanUpset.pdf", height = 5, width = 11)
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
      sets.bar.color = c("#ECD98B", "#AAAAC2",  "#03875C", "#9A4C43", "#4E3227")[c(2, 1, 2, 5, 4, 3, 4, 1, 5, 3)])
dev.off()

# human and sequin
DTU.gene.illumina.100vs000 <- DTU.gene.human.illumina.100vs000
DTU.gene.ONT.100vs000 <- DTU.gene.human.ONT.100vs000
for(i in 1:5){
  DTU.gene.illumina.100vs000[[i]] <- c(DTU.gene.illumina.100vs000[[i]], DTU.gene.sequin.illumina.100vs000[[i]])
  DTU.gene.ONT.100vs000[[i]] <- c(DTU.gene.ONT.100vs000[[i]], DTU.gene.sequin.ONT.100vs000[[i]])
}
pdf("plots/DTU/DTUgeneUpset.pdf", height = 5, width = 11)
upset(fromList(append(DTU.gene.illumina.100vs000, DTU.gene.ONT.100vs000)), 
      nsets=10, nintersects = 25, order.by = "freq",
      text.scale = c(1.5, 1.5, 1.5, 1.2, 1.2, 1.5),
      sets.bar.color = rep(c("#ECD98B", "#AAAAC2",  "#03875C", "#9A4C43", "#4E3227"), 2)[order(sapply(append(DTU.gene.illumina.100vs000, DTU.gene.ONT.100vs000), length, simplify=T), decreasing = T)])
dev.off()

# gene tested by all methods only
gene.human.ONT <- readRDS(file.path(DIR, "ONT/DTUmix/gene.human.RDS"))
gene.sequin.ONT <- readRDS(file.path(DIR, "ONT/DTUmix/gene.sequin.RDS"))
gene.human.illumina <- readRDS(file.path(DIR, "illumina/DTUmix/gene.human.RDS"))
gene.sequin.illumina <- readRDS(file.path(DIR, "illumina/DTUmix/gene.sequin.RDS"))
gene.human.ONT.100vs000 <- lapply(gene.human.ONT, function(x){
  return(x[[1]])
})
gene.sequin.ONT.100vs000 <- lapply(gene.sequin.ONT, function(x){
  return(x[[1]])
})
gene.human.illumina.100vs000 <- lapply(gene.human.illumina, function(x){
  return(x[[1]])
})
gene.sequin.illumina.100vs000 <- lapply(gene.sequin.illumina, function(x){
  return(x[[1]])
})
gene.human <- Reduce(intersect, append(gene.human.ONT.100vs000, gene.human.illumina.100vs000))
gene.sequin <- Reduce(intersect, append(gene.sequin.ONT.100vs000, gene.sequin.illumina.100vs000))
DTU.gene.illumina.100vs000.filt <- lapply(DTU.gene.illumina.100vs000, function(x){
  return(x[x %in% c(gene.human, gene.sequin)])
})
DTU.gene.ONT.100vs000.filt <- lapply(DTU.gene.ONT.100vs000, function(x){
  return(x[x %in% c(gene.human, gene.sequin)])
})
pdf("plots/DTU/DTUgeneUpsetFilt.pdf", height = 5, width = 11)
upset(fromList(append(DTU.gene.illumina.100vs000.filt, DTU.gene.ONT.100vs000.filt)), 
      nsets=10, nintersects = 25, order.by = "freq",
      text.scale = c(1.5, 1.5, 1.5, 1.2, 1.2, 1.5),
      sets.bar.color = rep(c("#ECD98B", "#AAAAC2",  "#03875C", "#9A4C43", "#4E3227"), 2)[order(sapply(append(DTU.gene.illumina.100vs000.filt, DTU.gene.ONT.100vs000.filt), length, simplify=T), decreasing = T)])
dev.off()

# some further exploration
DTU.tx.human.ONT.100vs000 <- lapply(DTU.tx.human.ONT, function(x){
  return(x[[1]])
})
DTU.tx.sequin.ONT.100vs000 <- lapply(DTU.tx.sequin.ONT, function(x){
  return(x[[1]])
})
DTU.tx.human.illumina.100vs000 <- lapply(DTU.tx.human.illumina, function(x){
  return(x[[1]])
})
DTU.tx.sequin.illumina.100vs000 <- lapply(DTU.tx.sequin.illumina, function(x){
  return(x[[1]])
})


DTU.ONT <- Reduce(union, append(DTU.tx.human.ONT.100vs000, DTU.tx.sequin.ONT.100vs000))
DTU.Illumina <- Reduce(union, append(DTU.tx.human.illumina.100vs000, DTU.tx.sequin.illumina.100vs000))
DTU.ONTonly <- DTU.ONT[!(DTU.ONT %in% DTU.Illumina)]
DTU.Illuminaonly <- DTU.Illumina[!(DTU.Illumina %in% DTU.ONT)]
DTU.ONTandIllumina <- intersect(DTU.Illumina, DTU.ONT)
library(edgeR)
s <- catchSalmon(file.path("../ONT/salmon_bs", list.files("../ONT/salmon_bs")))
dge <- DGEList(counts=s$counts/s$annotation$Overdispersion, genes=s$annotation)
s.short <- catchSalmon(file.path("../illumina/salmon_bs", list.files("../illumina/salmon_bs")))
dge.short <- DGEList(counts = s.short$counts/s.short$annotation$Overdispersion, genes = s.short$annotation)
rownames(dge) <- strsplit2(rownames(dge), "|", fixed = TRUE)[,1]
dge$genes$category <- "not DTU"
dge$genes$category[match(DTU.ONTonly, rownames(dge))] <- "ONT only"
dge$genes$category[match(DTU.Illuminaonly, rownames(dge))] <- "Illumina only"
dge$genes$category[match(DTU.ONTandIllumina, rownames(dge))] <- "ONT and Illumina"
rownames(dge.short) <- strsplit2(rownames(dge.short), "|", fixed = TRUE)[,1]
dge.short$genes$category <- "not DTU"
dge.short$genes$category[match(DTU.ONTonly, rownames(dge.short))] <- "ONT only"
dge.short$genes$category[match(DTU.Illuminaonly, rownames(dge.short))] <- "Illumina only"
dge.short$genes$category[match(DTU.ONTandIllumina, rownames(dge.short))] <- "ONT and Illumina"
dge$samples$group <- rep(c("000", "100", "075", "050", "025"), rep(3, 5))
dge <- dge[filterByExpr(dge),]
dge.short$samples$group <- rep(c("000", "100", "075", "050", "025"), rep(3, 5))
dge.short <- dge.short[filterByExpr(dge.short),]
dge$genes$totalCount <- rowSums(dge$counts)
dge.short$genes$totalCount <- rowSums(dge.short$counts)
dge$genes$category <- factor(dge$genes$category, levels = c("Illumina only", "ONT and Illumina", "ONT only", "not DTU"))
dge.short$genes$category <- factor(dge.short$genes$category, levels = c("Illumina only", "ONT and Illumina", "ONT only", "not DTU"))
pdf("plots/DTU/DTUcategory.pdf", height = 5, width = 8)
ggplot(dge$genes, aes(x=category, y=Length, fill=category)) +
  geom_violin() +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC", "gray80"))
ggplot(dge$genes, aes(x=category, y=EffectiveLength, fill=category)) +
  geom_violin() +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC", "gray80"))
ggplot(dge$genes, aes(x=category, y=totalCount, fill=category)) +
  geom_violin() +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC", "gray80"))+
  ggtitle("long read expression")
ggplot(dge.short$genes, aes(x=category, y=totalCount, fill=category)) +
  geom_violin() +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC", "gray80"))+
  ggtitle("short read expression")
ggplot(dge$genes, aes(x=category, y=Overdispersion, fill=category)) +
  geom_violin() +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC", "gray80"))+
  ggtitle("long read overdispersion")
ggplot(dge.short$genes, aes(x=category, y=Overdispersion, fill=category)) +
  geom_violin() +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC", "gray80"))+
  ggtitle("short read overdispersion")
dev.off()

# biotype and length of DTUs
txInfo.long <- readRDS("txInfo.long.RDS")
txInfo.short <- readRDS("txInfo.short.RDS")
DTU.gene.human.100vs000 <- data.frame(
  gene = unlist(append(DTU.gene.human.ONT.100vs000, DTU.gene.human.illumina.100vs000), use.names = TRUE),
  method = rep(rep(c("DEXSeq", "DRIMSeq", "edgeR", "limma", "satuRn"), 2),
               sapply(append(DTU.gene.human.ONT.100vs000, DTU.gene.human.illumina.100vs000), length, simplify = TRUE)),
  dataset = rep(rep(c("ONT", "Illumina"), c(5, 5)),
                sapply(append(DTU.gene.human.ONT.100vs000, DTU.gene.human.illumina.100vs000), length, simplify = TRUE))
)
library("GenomicFeatures")
gtf <- "/wehisan/home/allstaff/d/dong.x/annotation/HumanSequins/gencode.v33.sequins.gtf"
gtf_txdb <- makeTxDbFromGFF(gtf)
gene_list <- genes(gtf_txdb)
gene_list <- as.data.frame(gene_list)
saveRDS(gene_list, file="gene_list.RDS")
DTU.gene.human.100vs000$length <- gene_list[DTU.gene.human.100vs000$gene, "width"]
library(AnnotationHub)
ah <- AnnotationHub()
EnsDb.Hsapiens.v104 <- query(ah, c("EnsDb", "Homo Sapiens", 104))[[1]]
DTU.gene.human.100vs000$biotype <- mapIds(
  x = EnsDb.Hsapiens.v104,
  keys = substr(DTU.gene.human.100vs000$gene, 1, 15),
  keytype = "GENEID",
  column = "GENEBIOTYPE")

# deal with biotype
DTU.gene.human.100vs000$biotype[grepl("pseudogene$", DTU.gene.human.100vs000$biotype)] <- "pseudogene"
DTU.gene.human.100vs000 <- na.omit(DTU.gene.human.100vs000)

pdf("plots/DTU/DTUbiotypeGene.pdf", height = 5, width = 8)
ggplot(DTU.gene.human.100vs000, aes(x = method, fill=biotype))+
  geom_bar(position = "fill")+
  facet_grid(cols=vars(dataset)) +
  theme_bw() +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_fill_brewer(palette="Set2") +
  labs(fill = "Gene biotype", x = "Method", y = "Proportion of DTU gene")
dev.off()

library(ggridges)
pdf("plots/DTU/DTUlengthGene.pdf", height = 5, width = 8)
ggplot(DTU.gene.human.100vs000, aes(x = length, y=method, fill=method)) +
  geom_density_ridges(alpha = .7) +
  scale_fill_manual(values = c("#ECD98B", "#AAAAC2",  "#03875C", "#9A4C43", "#4E3227"))+
  scale_x_continuous(trans = "log10") +
  facet_grid(rows=vars(dataset)) +
  theme_bw() +
  theme(text = element_text(size = 20), legend.position = "NA")
dev.off()