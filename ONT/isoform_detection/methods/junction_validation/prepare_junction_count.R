# Using short-read junction count from SQANTI3 to validate long-read detected isoforms
# Prepare junction count data for plotting
# Author: Xueyi Dong

library(ggplot2)
library(ggridges)
library(ggsci)
library(parallel)

# prepare data----
load("ont_counts_filtered.RData")
junc_files <- list.files("./out", "junctions.txt", recursive = TRUE)
junc <- lapply(junc_files,
               function(x){
                 junc_count = read.delim(file.path("./out", x))
                 # method = strsplit2(x, "/")[1]
                 # filt = junc_count$isoform %in% eval(parse(text = 
                 #                                             paste0("ont_", method, "_counts")))[,1]
                 return(junc_count)
               })
names(junc) <- limma::strsplit2(junc_files, "/")[,1]
saveRDS(junc, "junc_unfiltered.RDS")
lapply(junc, function(x){
  length(unique(x$isoform))
})

# filter out lowly expressed transcripts
junc$bambu <- junc$bambu[junc$bambu$isoform %in% ont_bambu_counts$TXNAME, ]
flair_ids <- limma::strsplit2(ont_flair_counts$ids, "_")[ ,1]
junc$flair <- junc$flair[junc$flair$isoform %in% flair_ids, ]
flames_ids <- paste0("transcript:", ont_flames_counts$transcript_id)
junc$flames <- junc$flames[junc$flames$isoform %in% flames_ids, ]
junc$sqanti <- junc$sqanti[junc$sqanti$isoform %in% ont_sqanti_counts$Name, ]
junc$stringtie2 <- junc$stringtie2[junc$stringtie2$isoform %in% ont_stringtie_counts$Name, ]
junc$talon <- junc$talon[junc$talon$isoform %in% ont_talon_counts$annot_transcript_id, ]

saveRDS(junc, "junc.RDS")
