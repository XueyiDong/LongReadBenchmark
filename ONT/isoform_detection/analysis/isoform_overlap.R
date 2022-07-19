# Plot overlap of regions with common isoform coordinates
# Author: Mei Du

library(Biostrings)
library(rtracklayer)
library(IRanges)
library(UpSetR)
library(ComplexUpset)
library(ComplexHeatmap)
library(tidyverse)
library(here)

# load genome
gen <- readDNAStringSet(here("references","genome_rnasequin_decoychr_2.4.fa"))
gen_len <- width(gen)
names(gen_len) <- sub(" .*","",names(gen))
gen_bin <- Repitools::genomeBlocks(gen_len, width=2000) # make 2kb bins

# load annotations
annot_loc <- here("plots","readcount_cor","ONT_annotations")
bambu_annot <- import(paste0(annot_loc, "/ONT_extended_annotations_edited.gtf")) %>% .[.$type=="transcript",]
flair_annot <- import(paste0(annot_loc,"/flair.collapse.isoforms.gtf")) %>% .[.$type=="transcript",]
flames_annot <- import(paste0(annot_loc,"/isoform_annotated.filtered_flames.gff3")) %>% .[.$type=="transcript",]
seqlevelsStyle(flames_annot) <- "NCBI"
sqanti_annot <- import(paste0(annot_loc,"/merged_ont.collapsed_corrected.gtf")) %>% .[.$type=="transcript",]
stringtie_annot <- import(paste0(annot_loc,"/ONT_stringtie_merged_edited.gtf")) %>% .[.$type=="transcript",]
talon_annot <- import(paste0(annot_loc, "/ONT_filtered_annot_talon.gtf")) %>% .[.$type=="transcript",]

# load counts
load(here("plots","readcount_cor", "counts", "ont_counts_filtered.RData"))

# filter annotations
bambu_annot_filtered <- bambu_annot[bambu_annot$transcript_id %in% ont_bambu_counts$TXNAME,]
flair_annot_filtered <- flair_annot[flair_annot$transcript_id %in% str_extract(ont_flair_counts$ids, "[^_]*"),]
flames_annot_filtered <- flames_annot[flames_annot$transcript_id %in% ont_flames_counts$transcript_id,]
sqanti_annot_filtered <- sqanti_annot[sqanti_annot$transcript_id %in% ont_sqanti_counts$Name,]
stringtie_annot_filtered <- stringtie_annot[stringtie_annot$transcript_id %in% ont_stringtie_counts$Name,]
talon_annot_filtered <- talon_annot[talon_annot$transcript_id %in% ont_talon_counts$annot_transcript_id,]

# get genome-isoform overlaps
bambu_overlap <- queryHits(findOverlaps(gen_bin, bambu_annot_filtered))
flair_overlap <- queryHits(findOverlaps(gen_bin, flair_annot_filtered))
flames_overlap <- queryHits(findOverlaps(gen_bin, flames_annot_filtered))
sqanti_overlap <- queryHits(findOverlaps(gen_bin, sqanti_annot_filtered))
stringtie_overlap <- queryHits(findOverlaps(gen_bin, stringtie_annot_filtered))
talon_overlap <- queryHits(findOverlaps(gen_bin, talon_annot_filtered))

# binary matrix
overlap <- list(bambu = bambu_overlap,
                FLAIR = flair_overlap,
                FLAMES = flames_overlap,
                SQANTI3 = sqanti_overlap,
                StringTie2 = stringtie_overlap,
                TALON = talon_overlap)

# plot
tool <- c("bambu","FLAIR","FLAMES","SQANTI3","StringTie2","TALON")
bin_mat <- fromList(overlap)
upset <- ComplexUpset::upset(bin_mat[,order(ncol(bin_mat):1)], 
                             tool[order(length(tool):1)], 
                             n_intersections=20, 
                             width_ratio=0.2,
                             sort_sets=FALSE,
                             base_annotations = list(
                               'Intersection size'=(
                                 intersection_size()
                                 + theme(plot.background = element_blank(),
                                         panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank(),
                                         axis.title = element_text(size=13),
                                         axis.text = element_text(size=13)))
                               + ylab('Regions with common isoform coordinates')),
                             themes=upset_modify_themes(
                               list('intersections_matrix'=theme(text=element_text(size=16), axis.title.x = element_blank()),
                                    'overall_sizes'=theme(axis.text.x=element_text(size=15), axis.title.x = element_text(size=15)))))

#pdf(here("plots","isoform_overlap.pdf"),width=16, height=6)
upset
#dev.off()


