# Using short-read junction count from SQANTI3 to validate long-read detected isoforms
# Author: Xueyi Dong

library(ggplot2)
library(ggridges)
library(ggsci)
library(parallel)
library(plyr)
library(dplyr)

# plot junction count distribution ---- 
junc <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/code/analysis/junc.RDS")
lapply(junc, function(x){
  length(unique(x$isoform))
})
# deal with redundant junction information
junc_no_dup <- lapply(junc, function(x){
  x$junction = paste(x$chrom, x$genomic_start_coord, x$genomic_end_coord, sep = "_")
  x <- x[!duplicated(x$junction), ]
})
junc_df <- Reduce(rbind, junc_no_dup)
junc_df$Tool <- rep(c("bambu", "FLAIR", "FLAMES", "Cupcake", "StringTie2", "TALON"), 
                    sapply(junc_no_dup, nrow, simplify = TRUE))

# plot percentage of novel junctions supported by Illumina
# defind >10 counts as supported
junc_df$supported_by_Illumina <- junc_df$total_coverage_unique > 10
# pdf("NovelJuncIlluminaSupport2.pdf", height = 5, width = 5)
plot <- junc_df %>% 
  dplyr::filter(junction_category == "novel") %>% 
  dplyr::count(Tool = factor(Tool), supported_by_Illumina = factor(supported_by_Illumina)) %>% 
  plyr::ddply(., .(Tool), transform, pct = n / sum(n)) %>% 
  # mutate(pct = prop.table(n)) %>% 
  ggplot(aes(x = Tool, y = n, fill = supported_by_Illumina, label = scales::percent(pct))) +
  geom_bar(stat = "identity") +
  geom_text(position = position_stack(vjust = 0.7), size = 4) +
  scale_fill_manual(values = c("#A8D1D1", "#FD8A8A")) +
  labs(y = "Novel junctions", fill = "Supported by Illumina\n(counts > 10)") +
  theme_bw()+
  theme(text = element_text(size = 15), legend.position = "bottom", 
        axis.text.x = element_text(angle = 30, hjust = 1,size=15))
# dev.off()
