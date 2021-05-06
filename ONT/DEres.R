library(ggplot2)
library(UpSetR)

res.human <- list.files("./DEres", pattern="Human.tsv$", all.files = TRUE)
res.human <- sapply(res.human, 
                    function(x){read.table(file.path("./DEres", x), stringsAsFactors = FALSE, header=TRUE, sep="\t")},
                    simplify=FALSE)
res.human <-lapply(res.human,
                   function(x){
                     colnames(x)=c("logFC", "FDR")
                     return(na.omit(x))})
DEgenes.human <- lapply(res.human,
                        function(x){
                          rownames(x)[x$FDR < 0.05]
                        })
names(DEgenes.human) <- gsub("res", "", names(DEgenes.human))
names(DEgenes.human) <- gsub("Human.tsv", "", names(DEgenes.human))
upset(fromList(DEgenes.human), order.by = "degree")

#-------------------sequin----------------------------

res.sequin <- list.files("./DEres", pattern="Sequin.tsv$", all.files = TRUE)
res.sequin <- sapply(res.sequin, 
                    function(x){read.table(file.path("./DEres", x), stringsAsFactors = FALSE, header=TRUE, sep="\t")},
                    simplify=FALSE)
res.sequin <-lapply(res.sequin,
                   function(x){
                     colnames(x)=c("logFC", "FDR")
                     return(na.omit(x))})
DEgenes.sequin <- lapply(res.sequin,
                        function(x){
                          rownames(x)[x$FDR < 0.05]
                        })
