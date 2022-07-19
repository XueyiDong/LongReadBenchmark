# Sequins analysis
# Author: Mei Du

library(tidyverse)
library(UpSetR)
library(ggpubr)
library(ggforce)
library(RColorBrewer)
library(here)

# load sequins
s <- as.data.frame(read_tsv("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/Mike_seqin/annotations/rnasequin_isoforms_2.4.tsv"))

# load annotations
annot_loc <- here("plots","readcount_cor","ONT_annotations")
bambu_annot <- read.delim(here(annot_loc, "ONT_extended_annotations_edited.gtf"), header = FALSE) 
flair_annot <- read.delim(here(annot_loc,"flair.collapse.isoforms.gtf"), header = FALSE)
flames_annot <- read.delim(here(annot_loc,"isoform_annotated.filtered_flames.gff3"), header = FALSE)
sqanti_annot <- read.delim(here(annot_loc,"merged_ont.collapsed_classification.txt"))
stringtie_annot <- read.delim(here(annot_loc,"ONT_stringtie_merged_edited.gtf"), header = FALSE)
talon_annot <- read.delim(here(annot_loc, "ONT_filtered_annot_talon.gtf"), header = FALSE)

# isolate sequins-related transcript annotations (chrIS)
s_bambu <- bambu_annot[bambu_annot$V1=="chrIS" & bambu_annot$V3=="transcript",c(1:3,9)]
s_flair <- flair_annot[flair_annot$V1=="chrIS" & flair_annot$V3=="transcript",c(1:3,9)]
s_flames <- flames_annot[flames_annot$V1=="chrIS" & flames_annot$V3=="transcript",c(1:3,9)]
s_sqanti <- sqanti_annot[sqanti_annot$chrom=="chrIS",1:8]
s_stringtie <- stringtie_annot[stringtie_annot$V1=="chrIS" & stringtie_annot$V3=="transcript",c(1:3,9)]
s_talon <- talon_annot[talon_annot$V1=="chrIS" & talon_annot$V3=="transcript",c(1:3,9)]  %>% separate(.,V9,into=c("annot_gene_id","annot_transcript_id"), sep=";") 
s_talon$annot_gene_id <- str_replace(s_talon$annot_gene_id, "gene_id ", "")
s_talon$annot_transcript_id <- str_replace(s_talon$annot_transcript_id, " transcript_id ", "")

# add counts to annotations
load(here("plots","readcount_cor", "counts", "ont_counts_filtered.RData"))
s_bambu$count <- ont_bambu_counts[match(str_replace(separate(s_bambu, V9, into=c("gene", "transcript"), sep=";")$transcript, " transcript_id ", ""),ont_bambu_counts$TXNAME),9]
s_bambu <- na.omit(s_bambu)
s_flair$count<- ont_flair_counts[pmatch(str_replace(separate(s_flair, V9, into=c("gene", "transcript"), sep=";")$transcript, " transcript_id ", ""),ont_flair_counts$ids),8]
s_flair <- na.omit(s_flair)
s_flames$count <- ont_flames_counts[match(str_replace(separate(s_flames, V9, into=c("id", "transcript","gene","support"), sep=";")$transcript, "transcript_id=",""), ont_flames_counts$transcript_id),9]
s_flames <- na.omit(s_flames)
s_sqanti$count <- ont_sqanti_counts[match(s_sqanti$isoform, ont_sqanti_counts$Name),5]
s_sqanti <- s_sqanti[!is.na(s_sqanti$count),]
s_stringtie$count <- ont_stringtie_counts[match(str_replace(separate(s_stringtie, V9, into=c("gene", "transcript", "cov", "FPKM", "TPM"), sep=";")$transcript, " transcript_id ", ""), ont_stringtie_counts$Name),5]
s_stringtie <- na.omit(s_stringtie)
s_talon <- inner_join(s_talon, ont_talon_counts, by = c("annot_gene_id", "annot_transcript_id"))[,-c(15:20)]
s_talon <- na.omit(s_talon)

# sequins discovered (out of 160 total)
s_known_bambu <- s_bambu[str_detect(s_bambu$V9, paste(c(s$NAME),collapse="|")),]
s_known_flair <- s_flair[str_detect(s_flair$V9, paste(c(s$NAME),collapse="|")),]
s_known_flames <- s_flames[str_detect(s_flames$V9, paste0(c(s$NAME),";",collapse="|")),]
s_known_sqanti <- unique(s_sqanti[str_detect(s_sqanti$associated_transcript, paste(c(s$NAME),collapse="|")),"associated_transcript"])
s_known_stringtie <- s_stringtie[str_detect(s_stringtie$V9, paste(c(s$NAME),collapse="|")),]
s_known_talon <- s_talon[s_talon$V2=="Sequin",]

# isolate non-sequin-related isoforms
s_novel_bambu <- s_bambu[!str_detect(s_bambu$V9, paste(c(s$NAME),collapse="|")),]
s_novel_flair <- s_flair[!str_detect(s_flair$V9, paste(c(s$NAME),collapse="|")),]
s_novel_flames <- s_flames[!str_detect(s_flames$V9, paste0(c(s$NAME),";",collapse="|")),]
s_novel_sqanti <- s_sqanti[s_sqanti$associated_transcript=="novel",]
s_novel_stringtie <- s_stringtie[!str_detect(s_stringtie$V9, paste(c(s$NAME),collapse="|")),]
s_novel_talon <- s_talon[s_talon$V2=="TALON",]

# summary table
summary <- data.frame(tool = c("bambu", "FLAIR", "FLAMES","SQANTI3", "StringTie2", "TALON"), 
                          seq_related_annot = c(nrow(s_bambu), nrow(s_flair), nrow(s_flames),nrow(s_sqanti),nrow(s_stringtie), nrow(s_talon)),
                          known_seq_160 = c(nrow(s_known_bambu), nrow(s_known_flair), nrow(s_known_flames), length(s_known_sqanti),nrow(s_known_stringtie), nrow(s_known_talon)),
                          novel_annot = c(nrow(s_novel_bambu), nrow(s_novel_flair), nrow(s_novel_flames), nrow(s_novel_sqanti),nrow(s_novel_stringtie), nrow(s_novel_talon)))

# bar plot
colnames(summary)[c(3,4)] <- c("known sequins", "artefactual isoforms")
data <- pivot_longer(summary[,c(1,3,4)], c("known sequins", "artefactual isoforms"), names_to = "category",values_to="count")

seq_all <- ggbarplot(data, "tool", "count",
                           fill = "category", color = "category",
                           width =0.9,position = position_dodge(0.8),
                           palette = brewer.pal(5, "Set1")[c(2,3,4)],
                           xlab = "Method", ylab = "Count") + theme(legend.title = element_blank(), legend.position = "bottom") +
  geom_hline(yintercept=160,linetype = 2, size=1.2, col = brewer.pal(5, "Set1")[3]) +
  scale_y_continuous(breaks = c(0,80,160,seq(240,12000,480)),expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  facet_zoom(ylim=c(0,160)) +
  theme(legend.text=element_text(size=15),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=15),
        axis.text.y = element_text(size=13), axis.title = element_text(size=15)) 

#pdf(here("plots","sequins_detected_full.pdf"),width=9, height=5.5)
seq_all
#dev.off()

# sequin rediscovery analysis using downsampled dataset
remove_seq_tr <- read.delim(here("references","seqremoved.txt"), header = FALSE) # removed sequin transcripts
gtf <- read.delim(here("references","genome_rnasequin_decoychr_2.4_edited.gtf"), header = FALSE) # reference
gtf_removedseq <- data.frame(gtf[str_detect(gtf$V9, paste(c(str_trim(remove_seq_tr$V1)), collapse="|")),], row.names=NULL) # modified reference after sequins removed

# load annotations
sub_annot_loc <- here("plots","readcount_cor","ONT_annotations","ONT_sub")
bambu_sub_annot <- read.delim(here(sub_annot_loc, "ONT_sub_10_extended_annotations.gtf"), header = FALSE) 
flair_sub_annot <- read.delim(here(sub_annot_loc,"flair.collapse.isoforms.gtf"), header = FALSE)
flames_sub_annot <- read.delim(here(sub_annot_loc,"isoform_annotated.filtered.gff3"), header = FALSE)
sqanti_sub_annot <- read.delim(here(sub_annot_loc, "merged_ont_sub.collapsed_corrected.gtf"), header = FALSE)
sqanti_sub_annot <- read.delim(here(sub_annot_loc,"merged_ont_sub.collapsed_classification.txt"))
stringtie_sub_annot <- read.delim(here(sub_annot_loc,"ONT_stringtie_merged_sub.gtf"), header = FALSE)
talon_sub_annot <- read.delim(here(sub_annot_loc, "ONT_sub_filtered_annot_talon.gtf"), header = FALSE)

# load counts
sub_counts_loc <- here("plots","readcount_cor","counts","ONT","ONT_sub")
sub_bambu_counts <- read.delim(paste0(sub_counts_loc,"/bambu_ONT_sub_10_counts_transcript.txt")) %>% data.frame(., total_A=rowSums(.[3:5]), total_B=rowSums(.[6:8]))
sub_flair_counts <- read.delim(paste0(sub_counts_loc,"/flair_ont_sub_counts_matrix.tsv")) %>% data.frame(., total_A=rowSums(.[2:4]), total_B=rowSums(.[5:7]))
sub_flames_counts <- read_csv(paste0(sub_counts_loc,"/flames_ont_sub_transcript_count.csv.gz")) %>% data.frame(., total_A=rowSums(.[3:5]), total_B=rowSums(.[6:8]))
sub_sqanti_counts <- read.delim(paste0(sub_counts_loc,"/ont_sub_sqanti.sf"))
sub_stringtie_counts <- read.delim(paste0(sub_counts_loc,"/ont_sub_stringtie.sf"))
sub_talon_counts <- read.delim(paste0(sub_counts_loc,"/ONT_sub_talon_abundance_filtered.tsv")) %>% data.frame(., total_A=rowSums(.[12:14]), total_B=rowSums(.[15:17]))

# isolate sequins-related transcript annotations (chrIS)
s_bambu_sub <- bambu_sub_annot[bambu_sub_annot$V1=="chrIS" & bambu_sub_annot$V3=="transcript",c(1:5,9)]
s_flair_sub <- flair_sub_annot[flair_sub_annot$V1=="chrIS" & flair_sub_annot$V3=="transcript",c(1:5,9)]
s_flames_sub <- flames_sub_annot[flames_sub_annot$V1=="chrIS" & flames_sub_annot$V3=="transcript",c(1:5,9)]
s_sqanti_sub <- sqanti_sub_annot[sqanti_sub_annot$chrom=="chrIS",1:8]
s_sqanti_sub$isoform <- str_replace(s_sqanti_sub$isoform, "transcript_id ", "")
s_sqanti_sub <- inner_join(s_sqanti_sub, sqanti_class_sub)
s_stringtie_sub <- stringtie_sub_annot[stringtie_sub_annot$V1=="chrIS" & stringtie_sub_annot$V3=="transcript",c(1:5,9)]
s_talon_sub <- talon_sub_annot[talon_sub_annot$V1=="chrIS" & talon_sub_annot$V3=="transcript",c(4,5,9)] %>% separate(., V9, into=c("1", "2", "3", "4", "5", "6", "7", "8", "9"), sep=";")

# add counts
s_bambu_sub[,c("count_A","count_B")] <- sub_bambu_counts[match(str_replace(separate(s_bambu_sub, V9, into=c("gene", "transcript"), sep=";")$transcript, " transcript_id ", ""),sub_bambu_counts$TXNAME),9:10]
t <- str_replace(separate(s_flair_sub, V9, into=c("gene", "transcript"), sep=";")$transcript, " transcript_id ", "")
g <- str_replace(separate(s_flair_sub, V9, into=c("gene", "transcript"), sep=";")$gene, "gene_id ", "")
s_flair_sub[,c("count_A","count_B")]<- sub_flair_counts[match(paste0(t,"_",g),sub_flair_counts$ids),8:9]
s_flames_sub[,c("count_A","count_B")] <- sub_flames_counts[match(str_replace(separate(s_flames_sub, V9, into=c("id", "transcript","gene","support"), sep=";")$transcript, "transcript_id=",""), sub_flames_counts$transcript_id),9:10]
s_sqanti_sub$count <- sub_sqanti_counts[match(s_sqanti_sub$isoform, sub_sqanti_counts$Name),5]
s_stringtie_sub$count <- sub_stringtie_counts[match(str_replace(separate(s_stringtie_sub, V9, into=c("gene", "transcript", "cov", "FPKM", "TPM"), sep=";")$transcript, " transcript_id ", ""), sub_stringtie_counts$Name),5]
s_talon_sub[,c("count_A","count_B")] <- sub_talon_counts[match(str_replace(s_talon_sub[,4], " transcript_id ", ""), sub_talon_counts$annot_transcript_id),18:19]

# filter counts
s_bambu_sub <- s_bambu_sub[rowSums(s_bambu_sub[,7:8]) >= 10,]
s_flair_sub <- s_flair_sub[rowSums(s_flair_sub[,7:8]) >= 10,]
s_flames_sub <- s_flames_sub[rowSums(s_flames_sub[,7:8]) >= 10,]
s_sqanti_sub <- s_sqanti_sub[s_sqanti_sub$count >= 10,]
s_stringtie_sub <- s_stringtie_sub[s_stringtie_sub$count >= 10,]
s_talon_sub <- s_talon_sub[rowSums(s_talon_sub[,12:13]) >= 10,]

s_known_bambu_sub <- s_bambu_sub[str_detect(s_bambu_sub$V9, paste(c(s$NAME),collapse="|")),]
s_known_flair_sub <- s_flair_sub[str_detect(s_flair_sub$V9, paste(c(s$NAME),collapse="|")) ,]
s_known_flames_sub <- s_flames_sub[str_detect(s_flames_sub$V9, paste0(c(s$NAME),";",collapse="|")),]
s_known_sqanti_sub <- unique(s_sqanti_sub[str_detect(s_sqanti_sub$associated_transcript, paste(c(s$NAME),collapse="|")),"associated_transcript"])
s_known_stringtie_sub <- s_stringtie_sub[str_detect(s_stringtie_sub$V9, paste(c(s$NAME) ,collapse="|")),]
s_known_talon_sub <- s_talon_sub[str_detect(s_talon_sub[,4], paste(c(s$NAME), collapse="|")),]

summary_sub <- data.frame(tool = c("bambu", "FLAIR", "FLAMES","SQANTI3", "StringTie2", "TALON"), 
                              seq_related_annot = c(nrow(s_bambu_sub), nrow(s_flair_sub), nrow(s_flames_sub),nrow(s_sqanti_sub),nrow(s_stringtie_sub), nrow(s_talon_sub)),
                              known_seq_120 = c(nrow(s_known_bambu_sub), nrow(s_known_flair_sub), nrow(s_known_flames_sub), length(s_known_sqanti_sub),nrow(s_known_stringtie_sub), nrow(s_known_talon_sub)),
                              recovered_seq_40 = c(23,16,12,21,17,13))
summary_sub$fd <- summary_sub$seq_related_annot - (summary_sub$known_seq_120 + summary_sub$recovered_seq_40)

# bar plot
colnames(summary_sub)[c(3:5)] <- c("sequins kept", "sequins removed", "artefactual isoforms")
data_sub <- pivot_longer(summary_sub[,c(1,3:5)], c("sequins kept", "sequins removed", "artefactual isoforms"), names_to = "category",values_to="count")

seq_recovered <- ggbarplot(data_sub, "tool", "count",
                           fill = "category", color = "category",
                           width =0.8,position = position_dodge(0.8),
                           palette = brewer.pal(5, "Set1")[c(2,3,4)],
                           xlab = "Method", ylab = "Count") + theme(legend.title = element_blank(), legend.position = "bottom") +
  geom_hline(yintercept=40,linetype = 2, size=1.2, col = brewer.pal(5, "Set1")[4]) +
  geom_hline(yintercept=120,linetype = 2, size=1.2, col = brewer.pal(5, "Set1")[3]) +
  scale_y_continuous(breaks = c(0,40,80,120,seq(240,2800,240)), expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  facet_zoom(ylim=c(0,120)) +
  theme(legend.text=element_text(size=15),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=15),
        axis.text.y = element_text(size=13), axis.title = element_text(size=15))

#pdf(here("plots","sequins_recovered.pdf"),width=9, height=5.5
seq_recovered
#dev.off()
