library(ggplot2)
library(UpSetR)

# for easy change of dataset
OUT="./DEres_0_1"

res.human <- list.files(OUT, pattern="Human.tsv$", all.files = TRUE)
res.human <- sapply(res.human, 
                    function(x){read.table(file.path(OUT, x), stringsAsFactors = FALSE, header=TRUE, sep="\t")},
                    simplify=FALSE)
res.human <-lapply(res.human,
                   function(x){
                     colnames(x)=c("Gene", "logFC", "FDR")
                     return(na.omit(x))})
saveRDS(res.human, file = file.path(OUT, "res.human.RDS"))
DEgenes.human <- lapply(res.human,
                        function(x){
                          x$Gene[x$FDR < 0.05]
                        })
names(DEgenes.human) <- gsub("res", "", names(DEgenes.human))
names(DEgenes.human) <- gsub("Human.tsv", "", names(DEgenes.human))
# Human De genes upset plot
pdf(file.path(OUT, "plots/humanDEupset.pdf"), height = 5)
upset(fromList(DEgenes.human), order.by = "freq")
dev.off()

#-------------------sequin----------------------------

res.sequin <- list.files(OUT, pattern="Sequin.tsv$", all.files = TRUE)
res.sequin <- sapply(res.sequin, 
                    function(x){read.table(file.path(OUT, x), stringsAsFactors = FALSE, header=TRUE, sep="\t")},
                    simplify=FALSE)
res.sequin <-lapply(res.sequin,
                   function(x){
                     colnames(x)=c("Gene", "logFC", "FDR")
                     return(na.omit(x))})
saveRDS(res.sequin, file=file.path(OUT, "res.sequin.RDS"))
DEgenes.sequin <- lapply(res.sequin,
                        function(x){
                          x$Gene[x$FDR < 0.05]
                        })
# Sequin DE genes upset plot
upset(fromList(DEgenes.sequin), order.by = "freq")
# Compare to annotation, calculate FDR and TPR
anno <- read.table("/wehisan/home/allstaff/d/dong.x/annotation/sequins/rnasequin_genes_2.4.tsv", header = TRUE, stringsAsFactors = FALSE)
anno$logFC <- log(anno$MIX_B / anno$MIX_A)
res.sequin <- data.frame(
  gene =  c(res.sequin[[1]]$Gene, res.sequin[[2]]$Gene, res.sequin[[3]]$Gene, res.sequin[[4]]$Gene, res.sequin[[5]]$Gene),
  logFC = c(res.sequin[[1]]$logFC, res.sequin[[2]]$logFC, res.sequin[[3]]$logFC, res.sequin[[4]]$logFC, res.sequin[[5]]$logFC),
  FDR = c(res.sequin[[1]]$FDR, res.sequin[[2]]$FDR, res.sequin[[3]]$FDR, res.sequin[[4]]$FDR, res.sequin[[5]]$FDR),
  method = rep(c("DESeq2", "EBSeq", "edgeR", "limma", "NOISeq"), c(nrow(res.sequin[[1]]), nrow(res.sequin[[2]]), nrow(res.sequin[[3]]), nrow(res.sequin[[4]]), nrow(res.sequin[[5]])))
)
res.sequin$logFC[res.sequin$method %in% c("EBSeq", "NOISeq")] <- -res.sequin$logFC[res.sequin$method %in% c("EBSeq", "NOISeq")] 
res.sequin$logFC_expected <- anno$logFC[match(res.sequin$gene, anno$NAME)]
# sequin logFC expected vs estimated
pdf(file.path(OUT, "plots/sequinlogFC.pdf"), height = 4)
ggplot(res.sequin, aes(x=logFC_expected, y=logFC, colour=method))+
  geom_point() +
  geom_smooth(alpha=0.5) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw()
dev.off()

# FDR = false positive / predicted positive
FDR <- as.data.frame(t(sapply(unique(res.sequin$method), function(x){
  fdr = nrow(res.sequin[res.sequin$method == x & res.sequin$FDR < 0.05 & res.sequin$logFC_expected == 0,]) / nrow(res.sequin[res.sequin$method == x & res.sequin$FDR < 0.05,] )
  return(c(x, fdr))
}, simplify=TRUE)))
colnames(FDR) <- c("Method", "False discovery rate")
FDR$`False discovery rate` <- as.numeric(FDR$`False discovery rate`)
pdf(file.path(OUT, "plots/FDR.pdf"), height = 4)
ggplot(FDR, aes(x=Method, y=`False discovery rate`)) +
  geom_bar(stat="identity") +
  theme_bw()
dev.off()
# TPR = true positive / positive
TPR <- as.data.frame(t(sapply(unique(res.sequin$method), function(x){
  tpr = nrow(res.sequin[res.sequin$method == x & res.sequin$FDR < 0.05 & res.sequin$logFC_expected != 0,]) / nrow(res.sequin[res.sequin$method == x & res.sequin$logFC_expected != 0,] )
  return(c(x, tpr))
}, simplify=TRUE)))
colnames(TPR) <- c("Method", "True positive rate")
TPR$`True positive rate` <- as.numeric(TPR$`True positive rate`)
pdf(file.path(OUT, "plots/TPR.pdf"), height = 4)
ggplot(TPR, aes(x=Method, y=`True positive rate`)) +
  geom_bar(stat="identity") +
  theme_bw()
dev.off()

