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
