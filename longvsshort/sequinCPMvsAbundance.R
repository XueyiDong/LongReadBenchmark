library(ggplot2)
library(cowplot)
p.long <- readRDS("../ONT/sequinCPMvsAbundancePlot.RDS")
p.short <- readRDS("../illumina/sequinCPMvsAbundancePlot.RDS")
leg <- get_legend(p.long + theme(legend.position = "bottom"))
p.long <- p.long + theme(legend.position = "none") +
  ggtitle("ONT sequins")
p.short <- p.short + theme(legend.position = "none") +
  ggtitle("Illumina sequins")
pdf("plots/sequinCPMvsAbundance.pdf", height = 5, width = 11)
plot_grid(p.long,
          p.short,
          leg,
          rel_heights = c(1, 0.2))
dev.off()