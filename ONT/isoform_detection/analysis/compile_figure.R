# Compile isoform detection and quantification figure
# Author: Mei Du 

library(cowplot)
library(here)

# load plot scripts
scripts <- list.files(path=here("code", "analysis"), full.names = TRUE) 

# build grid
## sequins
source(scripts[7])
top_row <- plot_grid(seq_all, seq_recovered, nrow=1, labels=c("a","b"),  label_size = 24, label_fontface = "plain", rel_widths=c(1,1.1), scale=0.92)

## isoform classification
source(scripts[2])
mid_1_row <- plot_grid(isoform_cc, labels = "c", label_size = 24, label_fontface = "plain", scale = 0.92)

## full vs downsampled isoform comparison and short-read junctions
source(scripts[3])
source(scripts[6])
mid_2_row <- plot_grid(iso_compare, no_fsm_2, labels = c("d","e"), label_size = 24, label_fontface = "plain", rel_widths = c(2.5, 1.4), scale = 0.92)

## isoform overlap
source(scripts[4])
bottom_row <- plot_grid(upset, labels = "f", label_size = 24, label_fontface = "plain", scale=0.92)

# compile figure
pdf(here("plots","isoform.pdf"), width=16.5, height=22)
plot_grid(top_row, mid_1_row, mid_2_row, bottom_row, ncol=1, rel_heights=c(1,0.85,0.95,0.95))
dev.off()
