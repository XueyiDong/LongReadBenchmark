library(ggplot2)
library(edgeR)
library(MetBrewer)
library(scales)

DIR="/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark"
# load DGE lists
# s <- catchSalmon(file.path(DIR, "ONT/salmon_bs", list.files(file.path(DIR, "/ONT/salmon_bs"))))
# dge <- DGEList(counts=s$counts/s$annotation$Overdispersion, genes=s$annotation)
# s.short <- catchSalmon(file.path(DIR, "illumina/salmon_bs", list.files(file.path(DIR, "illumina/salmon_bs"))))
# dge.short <- DGEList(counts = s.short$counts/s.short$annotation$Overdispersion, genes = s.short$annotation)

dge <- readRDS("dge.rds")
dge.short <- readRDS("dge.short.rds")

# organize read num stat 
read.stat <- data.frame(
  sample = rep(c(paste("H1975", 1:3, sep = "-"),
              paste("HCC827", c(1, 2, 5), sep = "-")), 2),
  raw_reads = c(46809934, 40748569, 42074278, 51044877, 48219467, 39145316,
               31073747, 34059318, 27974052, 31329213, 134000622, 30331521),
  mapped_reads = c(41054341, 35939259, 36437436, 46878103, 40629774, 36962945,
                   28999906, 31495290, 26126733, 29046678, 125097661, 26408681),
  read_counts = c(dge$samples$lib.size[1:6], dge.short$samples$lib.size[1:6]),
  dataset = rep(c("ONT", "Illumina"), c(6,6))
)

read.stat <- data.table::melt(read.stat, id.vars = c("sample", "dataset"))
# read num plot
pdf("plots/readNum.pdf", height = 5, width = 8)
ggplot(read.stat, aes(x=variable, y=value, fill=sample, label = value))+
  geom_bar(stat="identity") +
  facet_grid(cols=vars(dataset)) +
  geom_text(aes(label = label_number_si()(value)), size = 4, colour = "white", position = position_stack(vjust = 0.5)) +
  ylab("Number") +
  xlab("Category") +
  theme_bw() +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6))+
  scale_fill_manual(values = met.brewer("Troy", 6))
dev.off()

# gene biotype info
dge.human <- dge[grep("^ENST", rownames(dge)), ]
library(Homo.sapiens)
geneid <- strsplit2(rownames(dge.human$counts), "\\|")[, 2]
geneid <- strsplit2(geneid, "\\.")[,1]
genes <- select(Homo.sapiens, keys=geneid, columns=c("SYMBOL", "TXCHROM"), 
                keytype="ENSEMBL")
dge.human$genes <- cbind(dge.human$genes, genes[match(geneid, genes$ENSEMBL),])

library(AnnotationHub)
ah <- AnnotationHub()
EnsDb.Hsapiens.v98 <- query(ah, c("EnsDb", "Homo Sapiens", 98))[[1]]
biotype<- mapIds(
  x = EnsDb.Hsapiens.v98,
  # NOTE: Need to remove gene version number prior to lookup.
  keys = na.omit(dge.human$genes$SYMBOL),
  keytype = "SYMBOL",
  column = "GENEBIOTYPE")
dge.human$genes$biotype <- biotype[match(dge.human$genes$SYMBOL, names(biotype))]
head(dge.human$genes)
# deal with biotype
# http://asia.ensembl.org/info/genome/genebuild/biotypes.html
# scaRNA is a kind of snoRNA

dge.human$genes$biotype[grepl("pseudogene$", dge.human$genes$biotype)] <- "pseudogene"
dge.human$genes$biotype[grepl("^TR", dge.human$genes$biotype)] <- "TR_gene"
dge.human$genes$biotype[grepl("^IG", dge.human$genes$biotype)] <- "IG_gene"
dge.human$genes$biotype[dge.human$genes$biotype %in% c("miRNA", "misc_RNA", 
                                                       "piRNA", "rRNA", "siRNA",
                                                       "snRNA", "snoRNA", "scaRNA",
                                                       "tRNA", "vaultRNA"
                                                       )] <- "ncRNA"


# x.human$genes$sumCount <- rowSums(x.human$counts)
# biotypeSum <- aggregate(x.human$gene$sumCount, by=list(x.human$genes$biotype), FUN=sum, simplify=TRUE)
# # pseudogene count prop
# sum(biotypeSum$x[grep("pseudogene$", biotypeSum$Group.1)]) / sum(biotypeSum$x)
# # for each sample
# sapply(1:6, function(x){
#   typesum = aggregate(x.human$counts[,x], by=list(x.human$genes$biotype), FUN=sum, simplify=TRUE)
#   sum(typesum$x[grep("pseudogene$", typesum$Group.1)]) / sum(typesum$x)
# }, simplify=TRUE)