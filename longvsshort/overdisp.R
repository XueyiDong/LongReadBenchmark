library(edgeR)
library(ggplot2)
library(viridis)

s <- catchSalmon(file.path("../ONT/salmon_bs", list.files("../ONT/salmon_bs")))
dge <- DGEList(counts=s$counts/s$annotation$Overdispersion, genes=s$annotation)

s.short <- catchSalmon(file.path("../illumina/salmon_bs", list.files("../illumina/salmon_bs")))
dge.short <- DGEList(counts = s.short$counts/s.short$annotation$Overdispersion, genes = s.short$annotation)

# long vs short scatter plot
m <- match(rownames(dge.short), rownames(dge))
overdisp <- data.frame(
  Overdispersion.long = dge$genes$Overdispersion[m],
  Overdispersion.short = dge.short$genes$Overdispersion
)
ggplot(overdisp, aes(x = Overdispersion.short, y = Overdispersion.long)) +
  geom_point()
cor(overdisp$Overdispersion.long, overdisp$Overdispersion.short, use = "complete.obs")
# not correlated

# calculate number of tx per gene
calcTxNum <- function(genes, isSequin){
  geneid = rep(NA, length(genes))
  genes.sequin = strsplit2(genes[isSequin], "_", fixed = TRUE)
  geneid[isSequin] = paste(genes.sequin[,1], genes.sequin[,2], sep = "_")
  genes.human = strsplit2(genes[!isSequin], "|", fixed = TRUE)
  geneid[!isSequin] = genes.human[,2]
  txcount = table(geneid)
  txnum = as.numeric(txcount[match(geneid, names(txcount))])
  txnum
}
dge$genes$nTranscript <- calcTxNum(rownames(dge), grepl("^R", rownames(dge)))
dge.short$genes$nTranscript <- calcTxNum(rownames(dge.short), grepl("^R", rownames(dge.short)))

overdisp2 <- data.frame(
  Overdispersion = c(dge$genes$Overdispersion, dge.short$genes$Overdispersion),
  Length = c(dge$genes$Length, dge.short$genes$Length),
  Data = rep(c("ONT long read", "Illumina short read"), c(nrow(dge), nrow(dge.short))),
  Gene = c(rownames(dge), rownames(dge.short)),
  AveExpr = c(rowSums(dge$counts), rowSums(dge.short$counts)),
  numberTranscript = c(dge$genes$nTranscript, dge.short$genes$nTranscript)
)

# violin
ggplot(overdisp2, aes(x=Data, y=Overdispersion)) +
  geom_boxplot() +
  theme_bw()

# calculate number of tx per gene
calcTxNum <- function(genes, isSequin){
  geneid = rep(NA, length(genes))
  genes.sequin = strsplit2(genes[isSequin], "_", fixed = TRUE)
  geneid[isSequin] = paste(genes.sequin[,1], genes.sequin[,2], sep = "_")
  genes.human = strsplit2(genes[!isSequin], "|", fixed = TRUE)
  geneid[!isSequin] = genes.human[,2]
  txcount = table(geneid)
  txnum = as.numeric(txcount[match(geneid, names(txcount))])
  txnum
}
dge$genes$nTranscript <- calcTxNum(rownames(dge), grepl("^R", rownames(dge)))
dge.short$genes$nTranscript <- calcTxNum(rownames(dge.short), grepl("^R", rownames(dge.short)))
overdisp2$numberTranscript <- c(dge$genes$nTranscript, dge.short$genes$nTranscript)
# scatter plot with colour
pdf("plots/overdisp.pdf", height = 5, width = 8)
overdisp$numberTranscripts <- dge.short$genes$nTranscript
ggplot(overdisp, aes(x = Overdispersion.short, y = Overdispersion.long, colour=numberTranscripts)) +
  geom_point() +
  theme_bw() +
  scale_colour_viridis_c(trans="log10", name = "Number of transcripts per gene")
dev.off()

#stratify by number of transcripts per gene
maxnum <- max(overdisp2$numberTranscript)
overdisp2$nTxGroup <- Hmisc::cut2(overdisp2$numberTranscript, cuts = c(1, 2, 6, 11, 21, 51, maxnum))
pdf("plots/overdispBox.pdf", height = 5, width = 8)
ggplot(overdisp2, aes(x=nTxGroup, y=Overdispersion, fill=Data, colour=Data)) +
  geom_boxplot(alpha=0.4) +
  # geom_jitter() +
  labs(x = "Number of transcripts per gene") +
  scale_y_continuous(trans = "log10") +
  theme_bw()+
  theme(text = element_text(size = 20)) +
  scale_fill_manual(values = c("#438DAC", "#FCB344")) +
  scale_colour_manual(values = c("#438DAC", "#FCB344"))
dev.off()

pdf("plots/overdispLength.pdf", height = 5, width = 8)
ggplot(overdisp2[overdisp2$numberTranscript==1,], aes(x=Length, y=Overdispersion)) +
  stat_binhex() +
  scale_x_continuous(trans = "log10") +
  theme_bw()
dev.off()
pdf("plots/overdispExp.pdf", height = 4)
  ggplot(overdisp2[overdisp2$numberTranscript==1,], aes(x=AveExpr, y=Overdispersion)) +
    stat_binhex() +
    scale_x_continuous(trans = "log10") +
    theme_bw()
dev.off()


# explore ONT or Illumina only findings
DE.human.illumina.100vs000 <- readRDS("DE.human.illumina.100vs000.RDS")
DE.human.ONT.100vs000 <- readRDS( "DE.human.ONT.100vs000.RDS")
DE.sequin.illumina.100vs000 <- readRDS("DE.human.illumina.100vs000.RDS")
DE.sequin.ONT.100vs000 <- readRDS( "DE.human.ONT.100vs000.RDS")

intersect.human.ONT <- Reduce(intersect, DE.human.ONT.100vs000)
intersect.human.Illumina <- Reduce(intersect, DE.human.illumina.100vs000)
union.human.ONT <- Reduce(union, DE.human.ONT.100vs000)
union.human.Illumina <- Reduce(union, DE.human.illumina.100vs000)
DE.human.ONTonly <- intersect.human.ONT[!(intersect.human.ONT %in% union.human.Illumina)]
DE.human.Illuminaonly <- intersect.human.Illumina[!(intersect.human.Illumina %in% union.human.ONT)]
intersect.human.all <- intersect(intersect.human.ONT, intersect.human.Illumina)  

intersect.sequin.ONT <- Reduce(intersect, DE.sequin.ONT.100vs000)
intersect.sequin.Illumina <- Reduce(intersect, DE.sequin.illumina.100vs000)
union.sequin.ONT <- Reduce(union, DE.sequin.ONT.100vs000)
union.sequin.Illumina <- Reduce(union, DE.sequin.illumina.100vs000)
DE.sequin.ONTonly <- intersect.sequin.ONT[!(intersect.sequin.ONT %in% union.sequin.Illumina)]
DE.sequin.Illuminaonly <- intersect.sequin.Illumina[!(intersect.sequin.Illumina %in% union.sequin.ONT)]
intersect.sequin.all <- intersect(intersect.sequin.ONT, intersect.sequin.Illumina)

rownames(dge) <- strsplit2(rownames(dge), "|", fixed = TRUE)[,1]
dge$genes$category <- NA
dge$genes$category[match(DE.human.ONTonly, rownames(dge))] <- "ONT_only"
dge$genes$category[match(DE.human.Illuminaonly, rownames(dge))] <- "Illumina_only"
dge$genes$category[match(intersect.human.all, rownames(dge))] <- "Intersect_all"
dge$genes$category[match(DE.sequin.ONTonly, rownames(dge))] <- "ONT_only"
dge$genes$category[match(DE.sequin.Illuminaonly, rownames(dge))] <- "Illumina_only"
dge$genes$category[match(intersect.sequin.all, rownames(dge))] <- "Intersect_all"
rownames(dge.short) <- strsplit2(rownames(dge.short), "|", fixed = TRUE)[,1]
dge.short$genes$category <- NA
dge.short$genes$category[match(DE.human.ONTonly, rownames(dge.short))] <- "ONT_only"
dge.short$genes$category[match(DE.human.Illuminaonly, rownames(dge.short))] <- "Illumina_only"
dge.short$genes$category[match(intersect.human.all, rownames(dge.short))] <- "Intersect_all"
dge.short$genes$category[match(DE.sequin.ONTonly, rownames(dge.short))] <- "ONT_only"
dge.short$genes$category[match(DE.sequin.Illuminaonly, rownames(dge.short))] <- "Illumina_only"
dge.short$genes$category[match(intersect.sequin.all, rownames(dge.short))] <- "Intersect_all"
dge$samples$group <- rep(c("000", "100", "075", "050", "025"), rep(3, 5))
dge <- dge[filterByExpr(dge),]
dge.short$samples$group <- rep(c("000", "100", "075", "050", "025"), rep(3, 5))
dge.short <- dge.short[filterByExpr(dge.short),]
dge$genes$totalCount <- rowSums(dge$counts)
dge.short$genes$totalCount <- rowSums(dge.short$counts)

pdf("plots/DTEcategory.pdf", height = 5, width = 8)
ggplot(dge$genes, aes(x=category, y=Length, fill=category)) +
  geom_violin() +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC"))
ggplot(dge$genes, aes(x=category, y=Length, fill=category)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC"))
ggplot(dge$genes, aes(x=category, y=EffectiveLength, fill=category)) +
  geom_violin() +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC"))
ggplot(dge$genes, aes(x=category, y=EffectiveLength, fill=category)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC"))
ggplot(dge$genes, aes(x=category, y=Overdispersion, fill=category)) +
  geom_violin() +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ggtitle("long read overdispersion") +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC"))
ggplot(dge.short$genes, aes(x=category, y=Overdispersion, fill=category)) +
  geom_violin() +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ggtitle("short read overdispersion") +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC"))
ggplot(dge$genes, aes(x=category, y=nTranscript, fill=category)) +
  geom_violin() +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC"))
ggplot(dge$genes, aes(x=category, y=totalCount, fill=category)) +
  geom_violin() +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ggtitle("long read expression") +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC"))
ggplot(dge$genes, aes(x=category, y=totalCount, fill=category)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ggtitle("long read expression") +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC"))
ggplot(dge.short$genes, aes(x=category, y=totalCount, fill=category)) +
  geom_violin() +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ggtitle("short read expression") +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC"))
ggplot(dge.short$genes, aes(x=category, y=totalCount, fill=category)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ggtitle("short read expression") +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC"))
dev.off()

dge.sequin <- dge[grepl("^R", rownames(dge)),]
dge.short.sequin <- dge.short[grepl("^R", rownames(dge.short)),]
pdf("plots/DTEcategory.sequin.pdf", height = 5, width = 8)
ggplot(dge.sequin$genes, aes(x=category, y=Length, fill=category)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC"))
ggplot(dge.sequin$genes, aes(x=category, y=EffectiveLength, fill=category)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC"))
ggplot(dge.sequin$genes, aes(x=category, y=Overdispersion, fill=category)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ggtitle("long read overdispersion") +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC"))
ggplot(dge.short.sequin$genes, aes(x=category, y=Overdispersion, fill=category)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ggtitle("short read overdispersion") +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC"))
ggplot(dge.sequin$genes, aes(x=category, y=nTranscript, fill=category)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC"))
ggplot(dge.sequin$genes, aes(x=category, y=totalCount, fill=category)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ggtitle("long read expression") +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC"))
ggplot(dge.short.sequin$genes, aes(x=category, y=totalCount, fill=category)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ggtitle("short read expression") +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC"))
dev.off()