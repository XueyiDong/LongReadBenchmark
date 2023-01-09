library(ggplot2)
library(MetBrewer)

qcdata <- readRDS("/vast/scratch/users/dong.x/long_read_benchmark/pilot_ONT/scripts/qc/summaryInfo.RDS")
qcdata <- as.data.frame(qcdata, stringsAsFactors = FALSE)
qcdata$Read_length <- as.numeric(qcdata$Read_length)
qcdata$Qscore <- as.numeric(qcdata$Qscore)
qcdata$Barcode[!(qcdata$Barcode %in% c(paste0("barcode0", 1:9), "barcode10"))] <- "other"

summary(qcdata[qcdata$Barcode != "other", "Read_length"])

# qcdata.sub <- qcdata[1:1000, ]
pdf("plots/LengthDistSample.pdf", height = 5, width = 8)
ggplot(qcdata, aes(x=Read_length, colour=Barcode, linetype=Barcode)) +
  geom_density() +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  scale_x_continuous(trans="log10") +
  labs(x= "ONT read length", colour="sample", linetype="sample") +
  scale_colour_manual(values = c(met.brewer("Troy", 10), "grey70")) +
  scale_linetype_manual(values = c(rep("solid", 10), "dashed"))
dev.off()

# undemultiplexed reads are filtered out
pdf("plots/LengthDist.pdf", height = 6, width = 5)
ggplot(qcdata[qcdata$Barcode != "other", ], aes(x=Read_length)) +
  geom_density(colour="black") +
  geom_histogram(aes(y=..density..), fill="#438DAC", colour="#438DAC", alpha = .3, bins=50) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  scale_x_continuous(trans="log10") +
  labs(x= "ONT read length")
dev.off()
