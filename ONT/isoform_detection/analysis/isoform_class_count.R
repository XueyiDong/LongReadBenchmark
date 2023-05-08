# Visualising SQANTI3-classified isoform structural categories and read counts 
# Authors: Xueyi Dong, Mei Du

library(tidyverse)
library(here)

# load SQANTI3 classifications
tool <- c("bambu", "flair", "flames", "sqanti", "stringtie", "talon")
for (i in 1:6) {
  assign(paste0("isoClass_",tool[i]),
         read.delim(list.files(path = "/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/ONT/isoform_detection/methods/junction_validation/out/", 
                    pattern = "classification.txt", 
                    recursive = TRUE,
                    full.names = TRUE)[i]))
}
isoClass_bambu$method <- "bambu"
isoClass_flair$method <- "FLAIR"
isoClass_flames$method <- "FLAMES"
isoClass_flames$isoform <- str_remove_all(isoClass_flames$isoform, "transcript:")
isoClass_sqanti$method <- "Cupcake"
isoClass_stringtie$method <- "StringTie2"
isoClass_talon$method <- "TALON"

# add counts
load(here("plots","readcount_cor", "counts", "ont_counts_filtered.RData"))

isoClass_bambu$count <- ont_bambu_counts[match(isoClass_bambu$isoform, ont_bambu_counts$TXNAME),"total"]
isoClass_flair$count <- ont_flair_counts[match(isoClass_flair$isoform, str_extract(ont_flair_counts$ids, "[^_]*")),"total"]
isoClass_flames$count <- ont_flames_counts[match(isoClass_flames$isoform, ont_flames_counts$transcript_id),"total"]
isoClass_sqanti$count <- ont_sqanti_counts[match(isoClass_sqanti$isoform, ont_sqanti_counts$Name),"NumReads"]
isoClass_stringtie$count <- ont_stringtie_counts[match(isoClass_stringtie$isoform, ont_stringtie_counts$Name),"NumReads"]
isoClass_talon$count <- ont_talon_counts[match(isoClass_talon$isoform, ont_talon_counts$annot_transcript_id),"total"]

isoClass <- plyr::rbind.fill(isoClass_bambu, isoClass_flair,isoClass_flames,isoClass_sqanti,isoClass_stringtie, isoClass_talon)
isoClass <- isoClass[!is.na(isoClass$count),]
isoClass$structural_category[isoClass$structural_category %in% c("genic", "genic_intron", "intergenic")] <- "genic/intergenic"
isoClass$structural_category <- factor(isoClass$structural_category, levels =c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "genic/intergenic", "fusion", "antisense"))

# plots
library(RColorBrewer)
library(scales)
library(ggplot2)
library(cowplot)
col.category <- brewer.pal(nlevels(isoClass$structural_category), "Set1")

plot_isoformClass <- ggplot(isoClass, aes(x = method, fill=structural_category)) + 
  geom_bar(position = position_stack(reverse = TRUE)) +
  theme_bw() +
  geom_hline(yintercept = 164, linetype="dashed") +
  labs(y = "Number of transcripts", x = "Method", fill = "Structural category") +
  theme(text = element_text(size = 20), legend.position = "none") +
  scale_fill_manual(values = col.category) +
  scale_y_continuous(labels = scales::comma)

my.labels <- c("bambu",
               "Cupcake+\nSalmon",
               "FLAIR", 
               "FLAMES",
               "StringTie2+\nSalmon",
               "TALON")


plot_isoformCount <- isoClass %>% dplyr::count(method, structural_category, wt = count, name = "count") %>% 
  ggplot() +
  geom_bar(aes(x=method, y=count, fill=structural_category), position = position_stack(reverse = TRUE), stat = "identity") +
  labs(y = "Read count", x = "Method", fill = "Structural category") +
  theme_bw() +
  theme(text = element_text(size = 20), legend.position = "bottom") +
  scale_fill_manual(values = col.category) +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6))+
  scale_x_discrete(labels= my.labels)

leg <- plot_grid(NULL,get_legend(plot_isoformCount), rel_widths = c(2.5,1))
isoform_cc <- plot_grid(plot_isoformClass, plot_isoformCount + theme(legend.position = "none"), leg, nrow=2,rel_heights = c(1, 0.2))

pdf(here("plots","isoform_ClassAndCount_new.pdf"), height = 5, width = 16)
isoform_cc
dev.off()
