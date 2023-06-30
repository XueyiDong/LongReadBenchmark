# A description based on distribution or histograms of the length, number of transcripts per gene, and number of exons would be more informative.
library(ggplot2)
library(cowplot)
library(dplyr)
library(plyranges)

# annodir <- "/wehisan/home/allstaff/d/dong.x/annotation/sequins"
annodir <- "/Volumes/dong.x/annotation/sequins"
anno <- read.delim(file.path(annodir, "rnasequin_isoforms_2.4.tsv"), sep = "\t", stringsAsFactors = FALSE)
geneanno <- read.delim(file.path(annodir, "rnasequin_genes_2.4.tsv"), sep = "\t", stringsAsFactors = FALSE)
gtf <- read_gff(file.path(annodir, "rnasequin_annotation_2.4.gtf"))

# gene length (including intron)
genes <- gtf %>% filter(type == "gene")
glen <- ggplot() +
  geom_histogram(aes(x = width(genes)), bins = 30)+
  scale_x_continuous(trans = "log10") +
  labs(x = "Gene length (including introns and exons)") +
  theme_bw()

# transcript length
tlen <- ggplot(anno, aes(x = LENGTH))+
  geom_histogram(bins = 30) +
  scale_x_continuous(trans = "log10") +
  labs(x = "Transcript length") +
  theme_bw()

# number of transcripts per gene
transcripts <- gtf %>% filter(type == "transcript")
txcounts <- as.data.frame(table(transcripts$gene_id))
tcount <- ggplot(txcounts, aes(x = Freq))+
  geom_bar() +
  labs(x = "Number of transcripts per gene") +
  theme_bw()

# number of exons per gene
exons <- gtf %>% filter(type == "exon")
exoncounts <- as.data.frame(table(exons$gene_id))
ecount <- ggplot(exoncounts, aes(x=Freq)) +
  geom_histogram(bins = 49)+
  labs(x = "Number of exons per gene") +
  theme_bw()

pdf("figures/supp/sequin.pdf", height = 5, width = 7)
plot_grid(glen, tlen, tcount, ecount, labels = c("a", "b", "c", "d"), ncol = 2)
dev.off()