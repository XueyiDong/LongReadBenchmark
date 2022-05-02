# Compile isoform detection and quantification figure
# Author: Mei Du 

library(cowplot)
library(here)

# load plot scripts
plot_dir <- here("plots")
scripts <- list.files(path=here("code", "analysis"),full.names = TRUE) 
scripts <- scripts[-1] # exclude this script

# build grid
source(scripts[5])
top_row <- plot_grid(seq_all, seq_recovered, nrow=1, labels=c("A","B"),  label_size = 18, rel_widths=c(1,1.1), scale=0.95)

source(scripts[3])
source(scripts[4])
mid_1_row <- plot_grid(hp, rp_filtered, nrow=1,rel_widths = c(2.5,1),labels = c("C","D"),label_size = 18,scale=0.95)

source(scripts[1])
mid_2_row <- plot_grid(isoform_cc, labels = "E", label_size = 18, scale=0.95)

source(scripts[2])
bottom_row <- plot_grid(upset,labels="F", label_size = 18, scale=0.95)

# compile figure
pdf(here(plot_dir,"/IsoformDetect_Quant.pdf"), width=16, height=22)
pdf("/stornext/General/data/user_managed/grpu_mritchie_0/Mei/long_read_benchmark/IsoformDetect_Quant.pdf", width=16.5, height=22)
plot_grid(top_row, mid_1_row, mid_2_row, bottom_row, ncol=1, rel_heights=c(1, 1.2, 0.85, 1))
dev.off()
