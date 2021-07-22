fc <- readRDS("counts_ONT.RDS")
x <- DGEList(counts = fc$counts)
x <- x[,1:6]
x$samples$group <- rep(c("H1975", "HCC827"), c(3,3))
x$genes <- fc$annotation
human.gene <- grep("^ENSG", rownames(x))
x.human <- x[human.gene, ]
library(Homo.sapiens)
geneid <- strsplit2(rownames(x.human$counts), "\\.")[,1]
genes <- select(Homo.sapiens, keys=geneid, columns=c("SYMBOL", "TXCHROM"), 
                keytype="ENSEMBL")
x.human$genes <- genes[match(geneid, genes$ENSEMBL),]
library(AnnotationHub)
ah <- AnnotationHub()
EnsDb.Hsapiens.v98 <- query(ah, c("EnsDb", "Homo Sapiens", 98))[[1]]
biotype<- mapIds(
  x = EnsDb.Hsapiens.v98,
  # NOTE: Need to remove gene version number prior to lookup.
  keys = na.omit(x.human$genes$SYMBOL),
  keytype = "SYMBOL",
  column = "GENEBIOTYPE")
x.human$genes$biotype <- biotype[match(x.human$genes$SYMBOL, names(biotype))]
head(x.human$genes)
x.human$genes$sumCount <- rowSums(x.human$counts)
biotypeSum <- aggregate(x.human$gene$sumCount, by=list(x.human$genes$biotype), FUN=sum, simplify=TRUE)
# pseudogene count prop
sum(biotypeSum$x[grep("pseudogene$", biotypeSum$Group.1)]) / sum(biotypeSum$x)
# for each sample
sapply(1:6, function(x){
  typesum = aggregate(x.human$counts[,x], by=list(x.human$genes$biotype), FUN=sum, simplify=TRUE)
  sum(typesum$x[grep("pseudogene$", typesum$Group.1)]) / sum(typesum$x)
}, simplify=TRUE)
