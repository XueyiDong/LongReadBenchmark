samples_ONT <- list.files("./ONT/salmon")
quant_ONT <- file.path("./ONT/salmon", samples_ONT, "quant.sf")
samples_Illumina <- list.files("./Illumina/salmon")
quant_Illumina <- file.path("./Illumina/salmon", samples_Illumina, "quant.sf")
library(tximport)
txi_ONT <- tximport(quant_ONT, type = "salmon", txOut=TRUE, countsFromAbundance = "no")
txi_Illumina <- tximport(quant_Illumina, type = "salmon", txOut=TRUE, countsFromAbundance = "no")
counts_ONT <- as.data.frame(txi_ONT$counts, stringAsFactors = FALSE)
counts_Illumina <- as.data.frame(txi_Illumina$counts, stringAsFactors = FALSE)
colnames(counts_ONT) <- samples_ONT
colnames(counts_Illumina) <- samples_Illumina
saveRDS(counts_ONT, file = "ONT/counts.RDS")
saveRDS(counts_Illumina, file = "Illumina/counts.RDS")
