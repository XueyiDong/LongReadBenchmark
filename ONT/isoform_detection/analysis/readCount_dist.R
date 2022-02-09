# Distribution of Illumina short-read TPM counts
# Author: Mei Du

library(ggplot2)
library(ggridges)
library(edgeR)
library(ggsci)
library(tidyverse)
library(here)

# load short-read counts
load(here("plots","readcount_cor", "counts", "ill_counts_filtered.RData"))

# log-TPM counts
ill_bambu_counts[,c("Illumina","tool")]<- data.frame(log2(ill_bambu_counts$TPM+0.01), "bambu")
ill_flair_counts[,c("Illumina","tool")]<- data.frame(log2(ill_flair_counts$TPM+0.01), "FLAIR")
ill_flames_counts[,c("Illumina","tool")]<- data.frame(log2(ill_flames_counts$TPM+0.01), "FLAMES")
ill_sqanti_counts[,c("Illumina","tool")]<- data.frame(log2(ill_sqanti_counts$TPM+0.01), "SQANTI3")
ill_stringtie_counts[,c("Illumina","tool")]<- data.frame(log2(ill_stringtie_counts$TPM+0.01), "StringTie2")
ill_talon_counts[,c("Illumina","tool")]<- data.frame(log2(ill_talon_counts$TPM+0.01), "TALON")

log_mat <- rbind(ill_bambu_counts, ill_flair_counts, ill_flames_counts, ill_sqanti_counts, ill_stringtie_counts, ill_talon_counts)
log_mat$tool <- as.factor(log_mat$tool)
log_mat_filtered <- log_mat[!log_mat$NumReads==0,] # remove 0 counts option

rp_filtered <- ggplot(log_mat_filtered, aes(x=Illumina, y = tool, fill = tool)) +
  geom_density_ridges(size = 0.1, scale=3.7) +
  guides(fill=guide_legend(reverse = TRUE)) +
  scale_y_discrete(expand=c(0,0))+
  scale_fill_jama(alpha = 0.3) +
  labs(x=expression("log"[2]*"-TPM Illumina read counts"), fill = "Tool") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_blank(),
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.line.x = element_line(size=0.1),axis.ticks.length=unit(0.25,"cm")) 
#pdf(here("plots","ont_illumina_readcount_dist_filtered_fig.pdf"),width=4, height=4)
rp_filtered
#dev.off()

