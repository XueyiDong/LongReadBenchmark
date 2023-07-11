library(edgeR)
library(ggplot2)
library(viridis)
library(dplyr)

DIR="/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark"

s <- catchSalmon(file.path(DIR, "ONT/salmon_bs", list.files(file.path(DIR, "/ONT/salmon_bs"))))
dge <- DGEList(counts=s$counts, genes=s$annotation)

s.short <- catchSalmon(file.path(DIR, "illumina/salmon_bs", list.files(file.path(DIR, "illumina/salmon_bs"))))
dge.short <- DGEList(counts = s.short$counts, genes = s.short$annotation)

# long vs short scatter plot ----
m <- match(rownames(dge.short), rownames(dge))
overdisp <- data.frame(
  Overdispersion.long = dge$genes$Overdispersion[m],
  Overdispersion.short = dge.short$genes$Overdispersion
)
ggplot(overdisp, aes(x = Overdispersion.short, y = Overdispersion.long)) +
  geom_point()

# calculate number of tx per gene ----
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
  Data = rep(c("ONT", "Illumina"), c(nrow(dge), nrow(dge.short))),
  Gene = c(rownames(dge), rownames(dge.short)),
  AveExpr = c(rowSums(dge$counts), rowSums(dge.short$counts)),
  numberTranscript = c(dge$genes$nTranscript, dge.short$genes$nTranscript)
)
overdisp2$numberTranscript <- c(dge$genes$nTranscript, dge.short$genes$nTranscript)

# # scatter plot with colour
# pdf("plots/overdisp.pdf", height = 5, width = 8)
# overdisp$numberTranscripts <- dge.short$genes$nTranscript
# ggplot(overdisp, aes(x = Overdispersion.short, y = Overdispersion.long, colour=numberTranscripts)) +
#   geom_point() +
#   theme_bw() +
#   scale_colour_viridis_c(trans="log10", name = "Number of transcripts per gene")
# dev.off()

#stratify by number of transcripts per gene
maxnum <- max(overdisp2$numberTranscript)
overdisp2$nTxGroup <- Hmisc::cut2(overdisp2$numberTranscript, cuts = c(1, 2, 6, 11, 21, 51, maxnum))
sample_sizes <- overdisp2 %>% group_by(nTxGroup, Data) %>% summarise(n = n(), max_y = max(Overdispersion))
# overdisp box plot ----
pdf("plots/overdispBoxNew.pdf", height = 4, width = 8)
p.ovd <- ggplot(overdisp2, aes(x=nTxGroup, y=Overdispersion, fill=Data, colour=Data)) +
  geom_boxplot(varwidth = TRUE, alpha=0.4, outlier.shape = NA) +
  geom_text(data = sample_sizes, aes(x = nTxGroup, y = max_y, colour = Data, label = paste0("n = ", n)),
            vjust = -0.5, hjust = ifelse(sample_sizes$Data == "ONT", -0.2, 1.2)) +
  # geom_violin(alpha=0) +
  # geom_jitter() +
  labs(x = "Number of transcripts per gene", y = "Assignment ambiguity") +
  scale_y_continuous(trans = "log10") +
  theme_bw()+
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom") +
  scale_fill_manual(values = c("#FCB344", "#438DAC")) +
  scale_colour_manual(values = c("#FCB344", "#438DAC"))
plot(p.ovd)
dev.off()

# overdispersion vs length ----
maxLength <- max(overdisp2$Length)
overdisp2$lengthGroup <- Hmisc::cut2(overdisp2$Length, cuts = c(0, 500, 1000, 2000, maxLength))
pdf("plots/overdispLength.pdf", height = 5, width = 8)
ggplot(overdisp2, aes(x=lengthGroup, y=Overdispersion, fill=Data, colour=Data)) +
  geom_boxplot(varwidth = TRUE, alpha=0.4, outlier.shape = NA) +
  # geom_violin(alpha=0) +
  # geom_jitter() +
  labs(x = "Annotated transcript length", y = "Assignment ambiguity") +
  scale_y_continuous(trans = "log10") +
  theme_bw()+
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom") +
  scale_fill_manual(values = c("#FCB344", "#438DAC")) +
  scale_colour_manual(values = c("#FCB344", "#438DAC"))
dev.off()

# sequin overdisp ----
# overdisp.sequin <- overdisp[grep("^R", rownames(dge.short)), ]
m = grepl("^R", rownames(dge))
m2 = grepl("^R", rownames(dge.short))
overdisp2.sequin <- data.frame(
  Overdispersion = c(dge$genes$Overdispersion[m], dge.short$genes$Overdispersion[m2]),
  Length = c(dge$genes$Length[m], dge.short$genes$Length[m2]),
  Data = rep(c("ONT", "Illumina"), c(sum(m), sum(m2))),
  Gene = c(rownames(dge)[m], rownames(dge.short)[m2]),
  AveExpr = c(rowSums(dge$counts)[m], rowSums(dge.short$counts)[m2]),
  numberTranscript = c(dge$genes$nTranscript[m], dge.short$genes$nTranscript[m2])
)
sample_sizes <- overdisp2.sequin %>% group_by(numberTranscript, Data) %>% summarise(n = n(), max_y = max(Overdispersion))
pdf("plots/overdispBoxSequinNew.pdf", height = 4, width = 4)
p.ovd.se <- ggplot(overdisp2.sequin, aes(x=as.character(numberTranscript), y=Overdispersion, fill=Data, colour=Data)) +
  geom_boxplot(varwidth = TRUE, alpha=0.4) +
  geom_text(data = sample_sizes, aes(x = numberTranscript, y = max_y, colour = Data, label = paste0("n = ", n)),
            vjust = -0.5, hjust = ifelse(sample_sizes$Data == "ONT", -0.2, 1.2)) +
  labs(x = "Number of transcripts per gene", y = "Assignment ambiguity") +
  scale_y_continuous(trans = "log10") +
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "bottom") +
  scale_fill_manual(values = c("#FCB344", "#438DAC")) +
  scale_colour_manual(values = c("#FCB344", "#438DAC"))
plot(p.ovd.se)
dev.off()

# human overdisp ----
overdisp2.human <- data.frame(
  Overdispersion = c(dge$genes$Overdispersion[!m], dge.short$genes$Overdispersion[!m2]),
  Length = c(dge$genes$Length[!m], dge.short$genes$Length[!m2]),
  Data = rep(c("ONT", "Illumina"), c(sum(!m), sum(!m2))),
  Gene = c(rownames(dge)[!m], rownames(dge.short)[!m2]),
  AveExpr = c(rowSums(dge$counts)[!m], rowSums(dge.short$counts)[!m2]),
  numberTranscript = c(dge$genes$nTranscript[!m], dge.short$genes$nTranscript[!m2])
)
maxnum <- max(overdisp2.human$numberTranscript)
overdisp2.human$nTxGroup <- Hmisc::cut2(overdisp2.human$numberTranscript, cuts = c(1, 2, 6, 11, 21, 51, maxnum))
p.ovd.hu <- ggplot(overdisp2.human, aes(x=nTxGroup, y=Overdispersion, fill=Data, colour=Data)) +
  geom_boxplot(varwidth = TRUE, alpha=0.4, outlier.shape = NA) +
  labs(x = "Number of transcripts per gene", y = "Assignment ambiguity") +
  scale_y_continuous(trans = "log10") +
  theme_bw()+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom") +
  scale_fill_manual(values = c("#FCB344", "#438DAC")) +
  scale_colour_manual(values = c("#FCB344", "#438DAC"))
plot(p.ovd.hu)

# plot overdispersion together ----
library(cowplot)
pdf("plots/overdispBoxAll2.pdf", height = 4, width = 8)
plot_grid(p.ovd.hu +
            ggtitle("Human transcripts") +
            geom_hline(yintercept = 10^layer_scales(p.ovd.se)$y$get_limits()[2],
                       alpha = 0.3,
                       linetype = "dashed") +
            geom_hline(yintercept = 1,
                       alpha = 0.3,
                       linetype = "dashed"),
          p.ovd.se +
            ggtitle("Sequins") +
            theme(axis.title.y = element_blank()),
          ncol = 2,
          rel_widths = c(2,1))
dev.off()


# explore ONT or Illumina only DE findings ----
DE.human.illumina.100vs000 <- readRDS("DE.human.illumina.100vs000.RDS")
DE.human.ONT.100vs000 <- readRDS( "DE.human.ONT.100vs000.RDS")
DE.sequin.illumina.100vs000 <- readRDS("DE.sequin.illumina.100vs000.RDS")
DE.sequin.ONT.100vs000 <- readRDS( "DE.sequin.ONT.100vs000.RDS")


intersect.human.ONT <- Reduce(intersect, DE.human.ONT.100vs000)
intersect.human.Illumina <- Reduce(intersect, DE.human.illumina.100vs000)
union.human.ONT <- Reduce(union, DE.human.ONT.100vs000)
union.human.Illumina <- Reduce(union, DE.human.illumina.100vs000)
DE.human.all <- intersect(union.human.ONT, union.human.Illumina)
DE.human.ONTonly <- union.human.ONT[!(union.human.ONT %in% union.human.Illumina)]
DE.human.Illuminaonly <- union.human.Illumina[!(union.human.Illumina %in% union.human.ONT)]

# # find DTE detected by all methods
# DE.human.all <- intersect(intersect.human.ONT, intersect.human.Illumina)
# # Find intersect DTE set by at least 2 method
# comb <- data.frame(
#   A=rep(1:5, rep(5, 5)),
#   B = rep(1:5, 5))
# comb <- comb[comb$A != comb$B,]
# union2.human.ONT <- apply(comb, 1, function(x){
#   return(intersect(DE.human.ONT.100vs000[[x[1]]], DE.human.ONT.100vs000[[x[2]]]))
# })
# union2.human.ONT <- Reduce(union, union2.human.ONT)
# DE.human.ONTonly <- union2.human.ONT[!(union2.human.ONT %in% union.human.Illumina)]
# union2.human.Illumina <- apply(comb, 1, function(x){
#   return(intersect(DE.human.illumina.100vs000[[x[1]]], DE.human.illumina.100vs000[[x[2]]]))
# })
# union2.human.Illumina <- Reduce(union, union2.human.Illumina)
# DE.human.Illuminaonly <- union2.human.Illumina[!(union2.human.Illumina %in% union.human.ONT)]

intersect.sequin.ONT <- Reduce(intersect, DE.sequin.ONT.100vs000)
intersect.sequin.Illumina <- Reduce(intersect, DE.sequin.illumina.100vs000)
union.sequin.ONT <- Reduce(union, DE.sequin.ONT.100vs000)
union.sequin.Illumina <- Reduce(union, DE.sequin.illumina.100vs000)
DE.sequin.all <- intersect(union.sequin.ONT, union.sequin.Illumina)
DE.sequin.ONTonly <- union.sequin.ONT[!(union.sequin.ONT %in% union.sequin.Illumina)]
DE.sequin.Illuminaonly <- union.sequin.Illumina[!(union.sequin.Illumina %in% union.sequin.ONT)]
# DE.sequin.all <- intersect(intersect.sequin.ONT, intersect.sequin.Illumina)
# union2.sequin.ONT <- apply(comb, 1, function(x){
#   return(intersect(DE.sequin.ONT.100vs000[[x[1]]], DE.sequin.ONT.100vs000[[x[2]]]))
# })
# union2.sequin.ONT <- Reduce(union, union2.sequin.ONT)
# DE.sequin.ONTonly <- union2.sequin.ONT[!(union2.sequin.ONT %in% union.sequin.Illumina)]
# union2.sequin.Illumina <- apply(comb, 1, function(x){
#   return(intersect(DE.sequin.illumina.100vs000[[x[1]]], DE.sequin.illumina.100vs000[[x[2]]]))
# })
# union2.sequin.Illumina <- Reduce(union, union2.sequin.Illumina)
# DE.sequin.Illuminaonly <- union2.sequin.Illumina[!(union2.sequin.Illumina %in% union.sequin.ONT)]

rownames(dge) <- strsplit2(rownames(dge), "|", fixed = TRUE)[,1]
dge$genes$category <- "Not DTE"
dge$genes$category[match(DE.human.ONTonly, rownames(dge))] <- "ONT only"
dge$genes$category[match(DE.human.Illuminaonly, rownames(dge))] <- "Illumina only"
dge$genes$category[match(DE.human.all, rownames(dge))] <- "ONT and Illumina"
dge$genes$category[match(DE.sequin.ONTonly, rownames(dge))] <- "ONT only"
dge$genes$category[match(DE.sequin.Illuminaonly, rownames(dge))] <- "Illumina only"
dge$genes$category[match(DE.sequin.all, rownames(dge))] <- "ONT and Illumina"
rownames(dge.short) <- strsplit2(rownames(dge.short), "|", fixed = TRUE)[,1]
dge.short$genes$category <- "Not DTE"
dge.short$genes$category[match(DE.human.ONTonly, rownames(dge.short))] <- "ONT only"
dge.short$genes$category[match(DE.human.Illuminaonly, rownames(dge.short))] <- "Illumina only"
dge.short$genes$category[match(DE.human.all, rownames(dge.short))] <- "ONT and Illumina"
dge.short$genes$category[match(DE.sequin.ONTonly, rownames(dge.short))] <- "ONT only"
dge.short$genes$category[match(DE.sequin.Illuminaonly, rownames(dge.short))] <- "Illumina only"
dge.short$genes$category[match(DE.sequin.all, rownames(dge.short))] <- "ONT and Illumina"
dge$samples$group <- rep(c("000", "100", "075", "050", "025"), rep(3, 5))
dge <- dge[filterByExpr(dge),]
dge.short$samples$group <- rep(c("000", "100", "075", "050", "025"), rep(3, 5))
dge.short <- dge.short[filterByExpr(dge.short),]
dge$genes$totalCount <- rowSums(dge$counts)
dge.short$genes$totalCount <- rowSums(dge.short$counts)
dge$genes$category <- factor(dge$genes$category, levels = c("Illumina only", "ONT and Illumina", "ONT only", "Not DTE"))
dge.short$genes$category <- factor(dge.short$genes$category, levels = c("Illumina only", "ONT and Illumina", "ONT only", "Not DTE"))
saveRDS(dge$genes, "DTEcategory.RDS")
saveRDS(dge.short$genes, "DTEcategory.short.RDS")

pdf("plots/DTE/DTEcategory.pdf", height = 5, width = 8)
ggplot(dge$genes, aes(x=category, y=Length, fill=category)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, outlier.colour = NA, alpha = 0) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC", "gray80"))+
  scale_colour_manual(values = c("#FCB344", "#A09F78", "#438DAC", "gray80"))

ggplot(dge$genes, aes(x=category, y=EffectiveLength, fill=category)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, outlier.colour = NA, alpha = 0) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC", "gray80"))

ggplot(dge$genes, aes(x=category, y=Overdispersion, fill=category)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, outlier.colour = NA, alpha = 0) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ggtitle("long read overdispersion") +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC", "gray80"))

ggplot(dge.short$genes, aes(x=category, y=Overdispersion, fill=category)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, outlier.colour = NA, alpha = 0) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ggtitle("short read overdispersion") +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC", "gray80"))

ggplot(dge$genes, aes(x=category, y=nTranscript, fill=category)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, outlier.colour = NA, alpha = 0) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC", "gray80"))

ggplot(dge$genes, aes(x=category, y=totalCount, fill=category)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, outlier.colour = NA, alpha = 0) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ggtitle("long read expression") +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC", "gray80"))

ggplot(dge.short$genes, aes(x=category, y=totalCount, fill=category)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, outlier.colour = NA, alpha = 0) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ggtitle("short read expression") +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC", "gray80"))
dev.off()

dge.sequin <- dge[grepl("^R", rownames(dge)),]
dge.short.sequin <- dge.short[grepl("^R", rownames(dge.short)),]
pdf("plots/DTE/DTEcategory.sequin.pdf", height = 5, width = 8)
ggplot(dge.sequin$genes, aes(x=category, y=Length, fill=category)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, outlier.colour = NA, alpha = 0) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC", "gray80"))+
  scale_colour_manual(values = c("#FCB344", "#A09F78", "#438DAC", "gray80"))

ggplot(dge.sequin$genes, aes(x=category, y=EffectiveLength, fill=category)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, outlier.colour = NA, alpha = 0) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC", "gray80"))

ggplot(dge.sequin$genes, aes(x=category, y=Overdispersion, fill=category)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, outlier.colour = NA, alpha = 0) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ggtitle("long read overdispersion") +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC", "gray80"))

ggplot(dge.short.sequin$genes, aes(x=category, y=Overdispersion, fill=category)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, outlier.colour = NA, alpha = 0) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ggtitle("short read overdispersion") +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC", "gray80"))

ggplot(dge.sequin$genes, aes(x=category, y=nTranscript, fill=category)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, outlier.colour = NA, alpha = 0) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC", "gray80"))

ggplot(dge.sequin$genes, aes(x=category, y=totalCount, fill=category)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, outlier.colour = NA, alpha = 0) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ggtitle("long read expression") +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC", "gray80"))

ggplot(dge.short.sequin$genes, aes(x=category, y=totalCount, fill=category)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, outlier.colour = NA, alpha = 0) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ggtitle("short read expression") +
  scale_fill_manual(values = c("#FCB344", "#A09F78", "#438DAC", "gray80"))
dev.off()