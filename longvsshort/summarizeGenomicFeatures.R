gene.counts.long <- readRDS("counts/gene.counts.long.RDS")
exon.counts.long <- readRDS("counts/exon.counts.long.RDS")
intron.counts.long <- readRDS("counts/intron.counts.long.RDS")
# intergenic.counts.long <- readRDS("counts/intergenic.counts.long.RDS")

gene.counts.short <- readRDS("counts/gene.counts.short.RDS")
exon.counts.short <- readRDS("counts/exon.counts.short.RDS")
intron.counts.short <- readRDS("counts/intron.counts.short.RDS")
# intergenic.counts.short <- readRDS("counts/intergenic.counts.short.RDS")

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
samples.long <- samples.short <- samples

# Calculate count number
samples.long$Total.Reads <- colSums(exon.counts.long$stat[,-1]) - unlist(exon.counts.long$stat[10, -1])
samples.long$Unmapped.reads <- unlist(exon.counts.long$stat[2, -1])
samples.long$Mapped <- samples.long$Total.Reads - samples.long$Unmapped.reads

# samples.long$subread.gene <- colSums(gene.counts.long$counts)
# samples.long$subread.exon <- colSums(exon.counts.long$counts)
# samples.long$subread.intron <- colSums(intron.counts.long$counts)
# samples.long$subread.intergenic <- colSums(intergenic.counts.long$counts)
samples.long$subread.gene <- unlist(gene.counts.long$stat[1, -1])
samples.long$subread.exon <- unlist(exon.counts.long$stat[1, -1])
samples.long$subread.intron <- samples.long$subread.gene - samples.long$subread.exon
samples.long$subread.intergenic <- samples.long$Mapped - samples.long$subread.gene

samples.short$Total.Reads <- colSums(exon.counts.short$stat[,-1])
samples.short$Unmapped.reads <- unlist(exon.counts.short$stat[2, -1])
samples.short$Mapped <- samples.short$Total.Reads - samples.short$Unmapped.reads

# samples.short$subread.gene <- colSums(gene.counts.short$counts)
# samples.short$subread.exon <- colSums(exon.counts.short$counts)
# samples.short$subread.intron <- colSums(intron.counts.short$counts)
# samples.short$subread.intergenic <- colSums(intergenic.counts.short$counts)
samples.short$subread.gene <- unlist(gene.counts.short$stat[1, -1])
samples.short$subread.exon <- unlist(exon.counts.short$stat[1, -1])
samples.short$subread.intron <- samples.short$subread.gene - samples.short$subread.exon
samples.short$subread.intergenic <- samples.short$Mapped - samples.short$subread.gene

samples.long <- data.table::melt(samples.long, id.vars = c("group", "sample"),
                                 variable.name = "Category")
samples.short <- data.table::melt(samples.short, id.vars = c("group", "sample"),
                                 variable.name = "Category")
samples.long$Dataset <- "ONT"
samples.short$Dataset <- "Illumina"
samples <- rbind(samples.long, samples.short)

write.csv(samples, file = "counts/features.mapping.summary.csv")

library(dplyr)
samples4plot <- samples %>% 
  filter(Category %in% c("Unmapped.reads", "subread.exon", "subread.intron", "subread.intergenic"))

library(ggplot2)
pdf("plots/featuredCountStat.pdf", height = 5, width = 8)
ggplot(samples4plot, aes(x=sample, y=value, fill=Category)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(cols = vars(Dataset)) +
  labs(y = "Proportion", x = "Sample") +
  theme_bw() +
  scale_fill_brewer(palette = "Dark2", 
                    labels = c("Unmapped reads", "Exonic reads", "Intronic reads", "Intergenic reads")) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()
