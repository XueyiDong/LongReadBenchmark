# A description based on distribution or histograms of the length, number of transcripts per gene, and number of exons would be more informative.
library(ggplot2)
library(cowplot)
library(dplyr)
library(plyranges)

annodir <- "/wehisan/home/allstaff/d/dong.x/annotation/sequins"
# annodir <- "/Volumes/dong.x/annotation/sequins"
anno <- read.delim(file.path(annodir, "rnasequin_isoforms_2.4.tsv"), sep = "\t", stringsAsFactors = FALSE)
geneanno <- read.delim(file.path(annodir, "rnasequin_genes_2.4.tsv"), sep = "\t", stringsAsFactors = FALSE)
gtf <- read_gff(file.path(annodir, "rnasequin_annotation_2.4.gtf"))

# gene length (including intron)
genes <- gtf %>% filter(type == "gene")
glen <- ggplot() +
  geom_histogram(aes(x = width(genes)), bins = 30)+
  scale_x_continuous(trans = "log10", limits = c(100, 3000000)) +
  labs(x = "Gene length (including introns and exons)") +
  theme_bw()
glen

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
  geom_histogram(binwidth = 1) +
  labs(x = "Number of transcripts per gene") +
  theme_bw()
tcount

# number of exons per gene
exons <- gtf %>% filter(type == "exon")
exons <- as.data.frame(exons)[, c("start", "end", "gene_id")]
exons <- exons[!duplicated(exons), ]
exoncounts <- as.data.frame(table(exons$gene_id))
ecount <- ggplot(exoncounts, aes(x=Freq)) +
  geom_histogram(binwidth = 1)+
  scale_x_continuous(limits = c(NA, 37))+
  labs(x = "Number of exons per gene") +
  theme_bw()
ecount

# expression fold change A vs B
anno$logFC <- log2(anno$MIX_A/anno$MIX_B)
lfc <- ggplot(anno, aes(x=logFC)) +
  geom_histogram() +
  labs(x = "log2 fold change of transcript expression") +
  theme_bw()
lfc

# usage FC
proportion <- readRDS("ONT/proportion.RDS")
propfc <- ggplot(proportion, aes(x=log2(FC)))+
  geom_histogram()+
  labs(x = "log2 fold change of transcript usage") +
  theme_bw()
propfc


pdf("figures/supp/sequin.pdf", height = 8, width = 8)
plot_grid(tcount, ecount, glen, tlen, lfc, propfc, labels = c("a", "b", "c", "d", "e", "f"), ncol = 2)
dev.off()

# plot_human <- readRDS("p_human.RDS")
# glen_human <- readRDS("p_human_glen.RDS")

plot_sequin <- plot_grid(tcount, ecount, glen,  ncol = 3, labels = c("a", "b", "c"))
saveRDS(plot_sequin, "plot_sequin.RDS")
plot_sequin
# plot_grid(plot_human, plot_sequin, ncol = 1, rel_heights = c(2, 1))