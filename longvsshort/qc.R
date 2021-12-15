library(ggplot2)
library(edgeR)

# load DGE lists
s <- catchSalmon(file.path("../ONT/salmon_bs", list.files("../ONT/salmon_bs")))
dge <- DGEList(counts=s$counts/s$annotation$Overdispersion, genes=s$annotation)
s.short <- catchSalmon(file.path("../illumina/salmon_bs", list.files("../illumina/salmon_bs")))
dge.short <- DGEList(counts = s.short$counts/s.short$annotation$Overdispersion, genes = s.short$annotation)

# organize read num stat 
read.stat <- data.frame(
  sample = rep(c(paste("H1975", 1:3, sep = "-"),
              paste("HCC827", c(1, 2, 5), sep = "-")), 2),
  raw_read = c(46809934, 40748569, 42074278, 51044877, 48219467, 39145316,
               31073747, 34059318, 27974052, 31329213, 134000622, 30331521),
  read_counts = c(dge$samples$lib.size[1:6], dge.short$samples$lib.size[1:6]),
  dataset = rep(c("ONT", "Illumina"), c(6,6))
)

read.stat <- data.table::melt(read.stat, id.vars = c("sample", "dataset"))
# read num plot
ggplot(read.stat, aes(x=sample, y=value, fill=variable))+
  geom_bar(stat="identity", position = "dodge") +
  facet_grid(cols=vars(dataset)) +
  theme_bw() +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1)) 

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
