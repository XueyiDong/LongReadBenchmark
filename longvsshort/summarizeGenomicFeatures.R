gene.counts.long <- readRDS("counts/gene.counts.long.RDS")
exon.counts.long <- readRDS("counts/exon.counts.long.RDS")
intron.counts.long <- readRDS("counts/intron.counts.long.RDS")
intergenic.counts.long <- readRDS("counts/intergenic.counts.long.RDS")

gene.counts.short <- readRDS("counts/gene.counts.short.RDS")
exon.counts.short <- readRDS("counts/exon.counts.short.RDS")
intron.counts.short <- readRDS("counts/intron.counts.short.RDS")
intergenic.counts.short <- readRDS("counts/intergenic.counts.short.RDS")

#Check that the file order is the same between all the objects

all(colnames(exon.counts.long$counts) == colnames(intron.counts.long$counts))
all(colnames(exon.counts.long$counts) == colnames(gene.counts.long$counts))
all(colnames(exon.counts.long$counts) == colnames(intergenic.counts.long$counts))

all(colnames(exon.counts.short$counts) == colnames(intron.counts.short$counts))
all(colnames(exon.counts.short$counts) == colnames(gene.counts.short$counts))
all(colnames(exon.counts.short$counts) == colnames(intergenic.counts.short$counts))

# Check total number of reads is the same
all(colSums(exon.counts.long$stat[,-1]) == colSums(intron.counts.long$stat[,-1]))
all(colSums(exon.counts.long$stat[,-1]) == colSums(gene.counts.long$stat[,-1]))
all(colSums(exon.counts.long$stat[,-1]) == colSums(intergenic.counts.long$stat[,-1]))

all(colSums(exon.counts.short$stat[,-1]) == colSums(intron.counts.short$stat[,-1]))
all(colSums(exon.counts.short$stat[,-1]) == colSums(gene.counts.short$stat[,-1]))
all(colSums(exon.counts.short$stat[,-1]) == colSums(intergenic.counts.short$stat[,-1]))

#Check number of unmapped reads is the same
all(exon.counts.long$stat[2, -1] == intron.counts.long$stat[2, -1])
all(exon.counts.long$stat[2, -1] == gene.counts.long$stat[2, -1])
all(exon.counts.long$stat[2, -1] == intergenic.counts.long$stat[2, -1])

all(exon.counts.short$stat[2, -1] == intron.counts.short$stat[2, -1])
all(exon.counts.short$stat[2, -1] == gene.counts.short$stat[2, -1])
all(exon.counts.short$stat[2, -1] == intergenic.counts.short$stat[2, -1])

#Check the order of samples
colnames(exon.counts.long$counts)
colnames(exon.counts.short$counts)

samples <- data.frame(
  group = rep(c("H1975", "HCC827"), each= 3),
  sample = c("H1975-1", "H1975-2", "H1975-3", "HCC827-1", "HCC827-2", "HCC827-5")
)

# Calculate count number
samples$Total.Reads.long <- colSums(exon.counts.long$stat[,-1])
samples$Unmapped.reads.long <- unlist(exon.counts.long$stat[2, -1])
samples$Mapped.long <- samples$Total.Reads.long - samples$Unmapped.reads.long

samples$subread.gene.long <- colSums(gene.counts.long$counts)
samples$subread.exon.long <- colSums(exon.counts.long$counts)
samples$subread.intron.long <- colSums(intron.counts.long$counts)
samples$subread.intergenic.long <- colSums(intergenic.counts.long$counts)

samples$Total.Reads.short <- colSums(exon.counts.short$stat[,-1])
samples$Unmapped.reads.short <- unlist(exon.counts.short$stat[2, -1])
samples$Mapped.short <- samples$Total.Reads.short - samples$Unmapped.reads.short

samples$subread.gene.short <- colSums(gene.counts.short$counts)
samples$subread.exon.short <- colSums(exon.counts.short$counts)
samples$subread.intron.short <- colSums(intron.counts.short$counts)
samples$subread.intergenic.short <- colSums(intergenic.counts.short$counts)

write.csv(samples, file = "counts/features.mapping.summary.csv")
