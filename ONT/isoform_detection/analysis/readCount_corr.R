# Long-read CPM and short-read TPM correlation
# Author: Mei Du

library(ggplot2)
library(edgeR)
library(tidyverse)
library(ggpubr)
library(here)

# load counts
counts_loc <- here("plots","readcount_cor","counts")

# ONT
ont_bambu_counts <- read.delim(paste0(counts_loc,"/ONT/bambu_ONT_counts_transcript.txt")) %>% data.frame(., total=rowSums(.[3:8]))
ont_flair_counts <- read.delim(paste0(counts_loc,"/ONT/flair_ont_counts_matrix.tsv")) %>% data.frame(., total=rowSums(.[2:7]))
ont_flames_counts <- read_csv(paste0(counts_loc,"/ONT/flames_ont_transcript_count.csv.gz")) %>% data.frame(., total=rowSums(.[3:8]))
ont_sqanti_counts <- read.delim(paste0(counts_loc,"/ONT/ont_sqanti.sf"))
ont_stringtie_counts <- read.delim(paste0(counts_loc,"/ONT/ont_stringtie.sf"))
ont_talon_counts <- read.delim(paste0(counts_loc,"/ONT/ONT_talon_abundance_filtered.tsv")) %>% data.frame(., total=rowSums(.[12:17]))

# Illumina
ill_bambu_counts <- read.delim(paste0(counts_loc,"/ONT/Illumina/bambu_ont_illumina.sf"))
ill_flair_counts <- read.delim(paste0(counts_loc,"/ONT/Illumina/flair_ont_illumina.sf"))
ill_flames_counts <- read.delim(paste0(counts_loc,"/ONT/Illumina/flames_ont_illumina.sf"))
ill_sqanti_counts <- read.delim(paste0(counts_loc,"/ONT/Illumina/sqanti_ont_illumina.sf"))
ill_stringtie_counts <- read.delim(paste0(counts_loc,"/ONT/Illumina/stringtie_ont_illumina.sf"))
ill_talon_counts <- read.delim(paste0(counts_loc,"/ONT/Illumina/talon_ont_illumina.sf"))

# filter transcripts
ont_bambu_counts <- ont_bambu_counts[ont_bambu_counts$total>=10,]
ont_flair_counts <- ont_flair_counts[ont_flair_counts$total>=10,]
ont_flames_counts <- ont_flames_counts[ont_flames_counts$total>=10,]
ont_sqanti_counts <- ont_sqanti_counts[ont_sqanti_counts$NumReads>=10,]
ont_stringtie_counts <- ont_stringtie_counts[ont_stringtie_counts$NumReads>=10,]
ont_talon_counts <- ont_talon_counts[ont_talon_counts$total>=10,]

# ensure short read counts match
ill_bambu_counts <- ill_bambu_counts[match(ont_bambu_counts$TXNAME, ill_bambu_counts$Name),]
ill_flair_counts <- ill_flair_counts[match(ont_flair_counts$ids, ill_flair_counts$Name),]
ill_flames_counts <- ill_flames_counts[match(ont_flames_counts$transcript_id, ill_flames_counts$Name),]
ill_sqanti_counts <- ill_sqanti_counts[match(ont_sqanti_counts$Name, ill_sqanti_counts$Name),]
ill_stringtie_counts <- ill_stringtie_counts[match(ont_stringtie_counts$Name, ill_stringtie_counts$Name),]
ill_talon_counts <- ill_talon_counts[match(ont_talon_counts$annot_transcript_id, ill_talon_counts$Name),]

# save counts
save(ont_bambu_counts, ont_flames_counts, ont_flair_counts, ont_sqanti_counts,ont_stringtie_counts,ont_talon_counts,
     file = paste0(counts_loc,"ont_counts_filtered.RData"))
save(ill_bambu_counts, ill_flames_counts, ill_flair_counts, ill_sqanti_counts,ill_stringtie_counts,ill_talon_counts,
     file = paste0(counts_loc,"ill_counts_filtered.RData"))

# read count correlation 
# count matrices
mat_bambu <- cbind(ill_bambu_counts, ont_bambu_counts)
mat_flair <- cbind(ill_flair_counts, ont_flair_counts)
mat_flames <- cbind(ill_flames_counts, ont_flames_counts)
mat_sqanti <- cbind(ill_sqanti_counts, ont_sqanti_counts)
colnames(mat_sqanti)[9:10] <- c("ontTPM","ontTotal")
mat_stringtie <- cbind(ill_stringtie_counts, ont_stringtie_counts)
colnames(mat_stringtie)[9:10] <- c("ontTPM","ontTotal")
mat_talon <- cbind(ill_talon_counts, ont_talon_counts) 

# log-transform matrices - short-read TPM, long-read CPM
log_mat_bambu <- data.frame(Illumina=log2(mat_bambu$TPM+0.5), ONT=cpm(mat_bambu$total, log=TRUE, prior.count=0.5), tool = "bambu")
log_mat_flair <- data.frame(Illumina=log2(mat_flair$TPM+0.5), ONT = cpm(mat_flair$total, log=TRUE, prior.count=0.5), tool = "FLAIR")
log_mat_flames <- data.frame(Illumina=log2(mat_flames$TPM+0.5), ONT = cpm(mat_flames$total, log=TRUE, prior.count=0.5), tool = "FLAMES")
log_mat_sqanti <- data.frame(Illumina=log2(mat_sqanti$TPM+0.5), ONT = cpm(mat_sqanti$ontTotal, log=TRUE, prior.count=0.5), tool = "SQANTI3")
log_mat_stringtie <- data.frame(Illumina=log2(mat_stringtie$TPM+0.5), ONT = cpm(mat_stringtie$ontTotal, log=TRUE, prior.count=0.5), tool = "StringTie2")
log_mat_talon <- data.frame(Illumina=log2(mat_talon$TPM+0.5), ONT = cpm(mat_talon$total, log=TRUE, prior.count=0.5), tool = "TALON")

log_mat <- rbind(log_mat_bambu, 
                 log_mat_flair, 
                 log_mat_flames, 
                 log_mat_sqanti,
                 log_mat_stringtie, 
                 log_mat_talon)
log_mat$tool <- as.factor(log_mat_ont$tool)

# hex bin scatter plot
# display feature number 
features <- data.frame(x=11.2, 
                       y= 16.5, 
                       label = c(paste("n =",nrow(log_mat_bambu)),
                                 paste("n =",nrow(log_mat_flair)),
                                 paste("n =",nrow(log_mat_flames)),
                                 paste("n =",nrow(log_mat_sqanti)),
                                 paste("n =",nrow(log_mat_stringtie)),
                                 paste("n =",nrow(log_mat_talon))),
                       tool = unique(log_mat_ont$tool))

# plot
hp <- ggplot(log_mat, aes(x=ONT,y=Illumina)) + 
  geom_hex(bins=50) +
  scale_fill_continuous(type="viridis", trans="log10") +
  xlab(expression("log"[2]*"CPM ONT read counts")) +
  ylab(expression("log"[2]*"TPM Illumina read counts")) +
  labs(fill=expression("log"[10]*" counts")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_blank(),axis.line = element_line(colour = "black")) +
  geom_text(data = features, aes(x,y,label=label)) +
  geom_smooth(method = "lm",se = FALSE,color = "orange") + stat_cor()+ border() + facet_wrap(.~tool)

#pdf(here("plots","readcount_cor","logCPM_hex_ont_illumina_fig.pdf"), height = 6, width = 10)
hp
#dev.off()

