# Compile isoform detection and quantification figure
# Author: Mei Du 

library(cowplot)
library(here)

# build grid
## sequins
source(here("code", "analysis","sequins_analysis.R"))
top_row <- plot_grid(seq_all,seq_sub, nrow=1, labels=c("a","b"),  label_size = 24, label_fontface = "bold", rel_widths=c(1.1,1), scale=0.92)

## isoform classification
source(here("code", "analysis","isoform_class_count.R"))
mid_1_row <- plot_grid(isoform_cc, labels = "c", label_size = 24, label_fontface = "bold", scale = 0.95)

## full vs downsampled isoform comparison and short-read junctions
source(c(here("code", "analysis","isoform_comparison.R"), here("code", "analysis","plot_junction_count.R")))
mid_2_row <- plot_grid(iso_compare, plot, labels = c("d","e"), label_size = 24, label_fontface = "bold", rel_widths = c(2.5, 1.4), scale = 0.92)

## isoform overlap
source(here("code", "analysis","isoform_overlap.R"))
bottom_row <- plot_grid(upset, labels = "f", label_size = 24, label_fontface = "bold", scale=0.92)

# compile figure
pdf("/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/plots/isoform_new3.pdf", width=16.5, height=23)
plot_grid(top_row, mid_1_row, mid_2_row, bottom_row, ncol=1, rel_heights=c(1.1,0.95,1,1))
dev.off()
