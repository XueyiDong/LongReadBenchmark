llibrary(ggplot2)
library(cowplot)
p.long <- readRDS("../ONT/sequinCPMvsAbundancePlot.RDS")
p.short <- readRDS("../illumina/sequinCPMvsAbundancePlot.RDS")
leg <- get_legend(p.long + theme(legend.position = "bottom")+
                    scale_colour_manual(name = "Sample", labels = paste(
                      rep(c("H1975", "HCC827"), c(3,3)),
                      c(1, 2, 3, 1, 2, 5),
                      sep = "-"
                    ), values = c("#E59EA1", "#DED8A2", "#F7C83E", "#4D5C28", "#DA4472", "#B77E60")))
p.long <- p.long + theme(legend.position = "none") +
  ggtitle("ONT sequins")
p.short <- p.short + theme(legend.position = "none") +
  ggtitle("Illumina sequins")
pdf("plots/sequinCPMvsAbundance.pdf", height = 8, width = 8)
plot_grid(p.long,
          p.short,
          leg,
          ncol = 1,
          rel_heights = c(1, 1, 0.2),
          byrow=FALSE)
dev.off()

pdf("plots/test.pdf")
p.long + theme(legend.position = "bottom") +
  scale_colour_manual(name = "Sample", labels = paste(
  rep(c("H1975", "HCC827"), c(3,3)),
  c(1, 2, 3, 1, 2, 5),
  sep = "-"
), values = c("#E59EA1", "#DED8A2", "#F7C83E", "#4D5C28", "#DA4472", "#B77E60"))
dev.off()