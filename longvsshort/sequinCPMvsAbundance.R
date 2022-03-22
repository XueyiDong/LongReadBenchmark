library(ggplot2)
library(cowplot)
library(MetBrewer)
p.long <- readRDS("../ONT/sequinCPMvsAbundancePlot.RDS")
p.short <- readRDS("../illumina/sequinCPMvsAbundancePlot.RDS")
leg <- get_legend(p.long + theme(legend.position = "right"))
p.long <- p.long + theme(legend.position = "none") +
  ggtitle("ONT sequins") 
p.short <- p.short + theme(legend.position = "none") +
  ggtitle("Illumina sequins") 
pdf("plots/sequinCPMvsAbundance.pdf", height = 8, width = 7)
plot_grid(p.long,
          p.short,
          leg,
          ncol = 2,
          # rel_heights = c(1, 1, 2),
          rel_widths = c(1, 0.3),
          byrow=FALSE)
dev.off()
