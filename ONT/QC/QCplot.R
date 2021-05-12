library(ggplot2)

# read number plot
load("../DE.RData")
stat <- as.data.frame(t(fc$stat))
colnames(stat) <- fc$stat[,1]
stat <- stat[-1,]
stat$sample <- rownames(stat)
stat$sample[1:6] <- paste0(rep(c("H1975_", "HCC827_"), c(3, 3)), rep(c("1", "2", "3"), 2))
pdf("plots/readNumber.pdf", height = 4)
ggplot(stat[c(1:6),], aes(x=sample, y=as.numeric(as.character(Assigned)))) + geom_bar(stat="identity") +
  labs(y="number of counts") +
  theme_bw() +  theme (axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# gene length vs expression plot
plotdata <- data.frame(
  Length=c(tt.human$Length, tt.sequin$Length),
  AveExpr = c(tt.human$AveExpr, tt.sequin$AveExpr),
  Source = rep(c("Human", "Sequin"),c(nrow(tt.human), nrow(tt.sequin)))
)
cor(plotdata$Length, plotdata$AveExpr)
pdf("plots/lengthCPM.pdf", height = 4)
ggplot(plotdata, aes(x=Length, y=AveExpr, colour=Source)) +
  geom_point(alpha=0.5) +
  scale_x_continuous(trans="log10") +
  labs(x="Gene length", y = "log2CPM") +
  theme_bw()
dev.off()
