library(ggplot2)

qcdata <- readRDS("./summaryInfo.RDS")
qcdata <- as.data.frame(qcdata, stringsAsFactors = FALSE)
qcdata$Read_length <- as.numeric(qcdata$Read_length)
qcdata$Qscore <- as.numeric(qcdata$Qscore)

maxLength = max(qcdata$Read_length)
maxLength
qcdata$LengthGroup <- Hmisc::cut2(qcdata$Read_length, cuts = c(0, 500, 1000,
                                                     3000, maxLength))
qcdata$LengthGroup <- gsub(" ", "", qcdata$LengthGroup)
qcdata$LengthGroup <- gsub(",", ", ", qcdata$LengthGroup)

qcdata$LengthGroup <- factor(qcdata$LengthGroup, levels = c(
  "[0, 500)", "[500, 1000)", "[1000, 3000)", "[3000, 208481]"
))

stat_box_data <- function(y, upper_limit = max(qcdata$Qscore) * 1.15) {
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = length(y)
    )
  )
}

pdf("plots/lengthQualityViolin.pdf", height = 4, width = 8)
ggplot(qcdata, aes(x=LengthGroup, y=Qscore, fill=LengthGroup, colour = LengthGroup)) +
  geom_violin(alpha = 0.4) +
  theme_bw() +
  theme(text = element_text(size=20), legend.position = "none") +
  labs(x = "Read length")+
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text", 
    hjust = 0.5,
    vjust = 0.9
  ) 
dev.off()

pdf("plots/sampleQualityBox.pdf", height = 4, width = 8)
ggplot(qcdata, aes(x=Barcode, y=Qscore, fill=Barcode, colour = Barcode)) +
  geom_boxplot(alpha = 0.4) +
  theme_bw() +
  theme(text = element_text(size=20), legend.position = "none", axis.text.x = element_text(angle = 90)) +
  labs(x = "Read length")+
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text", 
    hjust = 0.5,
    vjust = 0.9
  ) 
dev.off()
