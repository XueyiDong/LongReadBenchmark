library(ggplot2)
library(cowplot)
library(dplyr)
library(plyranges)


gtf <- "/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/HumanSequins/gencode.v33.sequins.gtf"
gtf <- read_gff(gtf)

# gene length (including intron)
genes <- gtf %>% filter(type == "gene")
glen <- ggplot() +
  geom_histogram(aes(x = width(genes)), bins = 30)+
  scale_x_continuous(trans = "log10", limits = c(100, 3000000)) +
  labs(x = "Gene length (including introns and exons)", title = "Human") +
  theme_bw()
glen

# number of transcripts per gene
transcripts <- gtf %>% filter(type == "transcript")
txcounts <- as.data.frame(table(transcripts$gene_id))
txcounts$highlight <- txcounts$Freq <= 4
tcount <- ggplot(txcounts, aes(x = Freq, fill = highlight))+
  geom_histogram(binwidth = 8) +
  labs(x = "Number of transcripts per gene", title = "Human") +
  scale_fill_manual(values = c("grey40", "indianred"))+
  theme_bw() +
  theme(legend.position = "none")
tcount

tcounts_hl <- ggplot(txcounts[txcounts$highlight, ], aes(x = Freq, fill=highlight))+
  geom_histogram(binwidth = 1)+
  labs(x = "Number of transcripts per gene", title = "Human") +
  scale_fill_manual(values = "indianred")+
  theme_bw() +
  theme(legend.position = "none")
tcounts_hl

# number of exons per gene
exons <- gtf %>% filter(type == "exon")
exons <- as.data.frame(exons)[, c("seqnames", "start", "end", "gene_id")]
exons <- exons[!duplicated(exons), ]
exoncounts <- as.data.frame(table(exons$gene_id))
exoncounts$highlight <- exoncounts$Freq <= 36
ecount <- ggplot(exoncounts, aes(x=Freq, fill=highlight)) +
  geom_histogram(binwidth = 8)+
  labs(x = "Number of exons per gene", title = "Human") +
  scale_fill_manual(values = c("grey40", "indianred"))+
  scale_x_continuous(limits = c(NA, 100)) +
  theme_bw() +
  theme(legend.position = "none")
ecount

ecount_hl <- ggplot(exoncounts[exoncounts$highlight, ], aes(x=Freq, fill=highlight)) +
  geom_histogram(binwidth = 1)+
  labs(x = "Number of exons per gene", title = "Human") +
  scale_fill_manual(values = "indianred")+
  scale_x_continuous(limits = c(NA, 37)) +
  theme_bw() +
  theme(legend.position = "none")
ecount_hl

# pdf("figures/supp/human.pdf", height = 3, width = 9)
# plot_grid(tcount, ecount, glen, labels = c("a", "b", "c"), ncol = 3)
# dev.off()

plot_human <- plot_grid(tcount, ecount, glen, tcounts_hl, ecount_hl,
                        labels = c("d", "e", "f", "g", "h"),
                        ncol = 3)
# saveRDS(plot_human, "p_human.RDS")
# saveRDS(glen, "p_human_glen.RDS")

plot_sequin <- readRDS("plot_sequin.RDS")
pdf("figures/supp/annoGeneInfo.pdf", height = 12, width = 12)
plot_grid(plot_sequin, plot_human, ncol = 1, rel_heights = c(1, 2))
dev.off()