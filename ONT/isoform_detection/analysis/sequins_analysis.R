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

########### full sequin recovery analysis ############
# annotations
ont_annot_loc <- "/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/plots/readcount_cor/ONT_annotations/"
bambu_ont_annot <- read.delim(paste0(ont_annot_loc, "ONT_extended_annotations_edited.gtf"), header = FALSE) 
flair_ont_annot <- read.delim(paste0(ont_annot_loc,"flair.collapse.isoforms.gtf"), header = FALSE)
flames_ont_annot <- read.delim(paste0(ont_annot_loc,"isoform_annotated.filtered_flames.gff3"), header = FALSE)
sqanti_ont_annot <- read.delim(paste0(ont_annot_loc,"merged_ont.collapsed_corrected.gtf"), header=FALSE) 
sqanti_class_ont <- read.delim(paste0(ont_annot_loc,"merged_ont.collapsed_classification.txt"))
stringtie_ont_annot <- read.delim(paste0(ont_annot_loc,"ONT_stringtie_merged_edited.gtf"), header = FALSE)
talon_ont_annot <- read.delim(paste0(ont_annot_loc, "ONT_filtered_annot_talon.gtf"), header = FALSE)

# isolate sequins-related transcript annotations - chrIS
s_bambu_ont <- bambu_ont_annot[bambu_ont_annot$V1=="chrIS" & bambu_ont_annot$V3=="transcript",c(1:3,9)]
s_flair_ont <- flair_ont_annot[flair_ont_annot$V1=="chrIS" & flair_ont_annot$V3=="transcript",c(1:3,9)]
s_flames_ont <- flames_ont_annot[flames_ont_annot$V1=="chrIS" & flames_ont_annot$V3=="transcript",c(1:3,9)]
s_sqanti_ont <- sqanti_ont_annot[sqanti_ont_annot$V1=="chrIS" & sqanti_ont_annot$V3=="transcript",] %>% separate(., V9, into=c("isoform", "gene"), sep = ";")
s_sqanti_ont$isoform <- str_replace(s_sqanti_ont$isoform, "transcript_id ", "")
s_sqanti_ont <- inner_join(s_sqanti_ont, sqanti_class_ont)
s_stringtie_ont <- stringtie_ont_annot[stringtie_ont_annot$V1=="chrIS" & stringtie_ont_annot$V3=="transcript",c(1:3,9)]
s_talon_ont <- talon_ont_annot[talon_ont_annot$V1=="chrIS" & talon_ont_annot$V3=="transcript",c(1:3,9)]  %>% separate(.,V9,into=c("annot_gene_id","annot_transcript_id"), sep=";") 
s_talon_ont$annot_gene_id <- str_replace(s_talon_ont$annot_gene_id, "gene_id ", "")
s_talon_ont$annot_transcript_id <- str_replace(s_talon_ont$annot_transcript_id, " transcript_id ", "")

# add counts
counts_loc <- "/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/plots/readcount_cor/counts/ONT/"
ont_bambu_counts <- read.delim(paste0(counts_loc,"bambu_ONT_counts_transcript.txt")) %>% data.frame(., total=rowSums(.[3:8]))
ont_flair_counts <- read.delim(paste0(counts_loc,"flair_ont_counts_matrix.tsv")) %>% data.frame(., total=rowSums(.[2:7]))
ont_flames_counts <- read_csv(paste0(counts_loc,"flames_ont_transcript_count.csv.gz")) %>% data.frame(., total=rowSums(.[3:8]))
ont_sqanti_counts <- read.delim(paste0(counts_loc,"/ont_sqanti.sf"))
ont_stringtie_counts <- read.delim(paste0(counts_loc,"ont_stringtie.sf"))
ont_talon_counts <- read.delim(paste0(counts_loc,"ONT_talon_abundance_filtered.tsv")) %>% data.frame(., total=rowSums(.[12:17]))

s_bambu_ont$count <- ont_bambu_counts[match(str_replace(separate(s_bambu_ont, V9, into=c("gene", "transcript"), sep=";")$transcript, " transcript_id ", ""),ont_bambu_counts$TXNAME),9]
s_flair_ont$count<- ont_flair_counts[pmatch(str_replace(separate(s_flair_ont, V9, into=c("gene", "transcript"), sep=";")$transcript, " transcript_id ", ""),ont_flair_counts$ids),8]
s_flair_ont <- na.omit(s_flair_ont)
s_flames_ont$count <- ont_flames_counts[match(str_replace(separate(s_flames_ont, V9, into=c("id", "transcript","gene","support"), sep=";")$transcript, "transcript_id=",""), ont_flames_counts$transcript_id),9]
s_sqanti_ont$count <- ont_sqanti_counts[match(s_sqanti_ont$isoform, ont_sqanti_counts$Name),5]
s_stringtie_ont$count <- ont_stringtie_counts[match(str_replace(separate(s_stringtie_ont, V9, into=c("gene", "transcript", "cov", "FPKM", "TPM"), sep=";")$transcript, " transcript_id ", ""), ont_stringtie_counts$Name),5]
s_talon_ont <- inner_join(s_talon_ont, ont_talon_counts, by = c("annot_gene_id", "annot_transcript_id"))[,-c(15:20)]

# filter counts
s_bambu_ont <- s_bambu_ont[s_bambu_ont$count >= 10,]
s_flair_ont <- s_flair_ont[s_flair_ont$count >= 10,]
s_flames_ont <- s_flames_ont[s_flames_ont$count >= 10,]
s_sqanti_ont <- s_sqanti_ont[s_sqanti_ont$count >= 10,]
s_stringtie_ont <- s_stringtie_ont[s_stringtie_ont$count >= 10,]
s_talon_ont <- s_talon_ont[s_talon_ont$total >= 10,]

# sequins discovered (out of 160)
s_known_bambu_ont <- s_bambu_ont[str_detect(s_bambu_ont$V9, paste(c(s$NAME),collapse="|")),]
s_known_flair_ont <- s_flair_ont[str_detect(s_flair_ont$V9, paste(c(s$NAME),collapse="|")),]
s_known_flames_ont <- s_flames_ont[str_detect(s_flames_ont$V9, paste0(c(s$NAME),";",collapse="|")),]
s_known_sqanti_ont <- unique(s_sqanti_ont[str_detect(s_sqanti_ont$associated_transcript, paste(c(s$NAME),collapse="|")),"associated_transcript"])
s_known_stringtie_ont <- s_stringtie_ont[str_detect(s_stringtie_ont$V9, paste(c(s$NAME),collapse="|")),]
s_known_talon_ont <- s_talon_ont[s_talon_ont$V2=="Sequin",]

# summary table
summary <- data.frame(tool = c("bambu", "Cupcake","FLAIR", "FLAMES","StringTie2","TALON"), 
                      seq_related_annot = c(nrow(s_bambu_ont), nrow(s_sqanti_ont),nrow(s_flair_ont), nrow(s_flames_ont),nrow(s_stringtie_ont), nrow(s_talon_ont)),
                      known_seq_160 = c(nrow(s_known_bambu_ont),length(s_known_sqanti_ont), nrow(s_known_flair_ont), nrow(s_known_flames_ont), nrow(s_known_stringtie_ont), nrow(s_known_talon_ont)))

summary$kept_fd <- summary$seq_related_annot - summary$known_seq_160
summary$precision <- summary$known_seq_160 / (summary$known_seq_160 + summary$kept_fd)
summary$recall <- summary$known_seq_160 / 160

# bar plot
data <- pivot_longer(summary, c("precision", "recall"), names_to = "rate",values_to="value")
seq_all <- ggbarplot(data, "tool", "value",
                     fill = "rate", color = "rate",
                     width =0.9,position = position_dodge(0.9),
                     palette = brewer.pal(5, "Set1")[c(2,3,4,5)],
                     xlab = "Method", ylab = "Rate") + 
  theme(axis.text.x = element_text( size=15),legend.title = element_blank(),legend.position="bottom",
        legend.text = element_text(size=15),
        axis.text.y = element_text(size=13), axis.title = element_text(size=15))

# pdf(here("sequins_detected_rates.pdf"),width=9, height=8)
seq_all
# dev.off()

############ sequin rediscovery analysis using downsampled dataset ############
# gtf_removedseq <- read_table(here("seq_removed.txt"),col_names=FALSE)
removed_genes <- c("R1_92", "R2_24", "R2_6", "R2_76", "R2_73", "R2_71", "R2_7", "R2_68",
                   "R2_47", "R2_55", "R2_46", "R2_41", "R2_42", "R2_38", "R2_32", "R2_26", "R2_153",
                   "R2_151", "R1_23", "R2_14", "R1_103", "R2_19", "R1_102", "R1_13","R1_101",
                   "R2_115", "R1_83", "R1_81", "R1_73", "R1_61", "R1_62", "R1_71", "R1_43", "R1_63",
                   "R1_51", "R1_41", "R1_42", "R1_11", "R2_57","R2_59")
kept_genes <- unique(str_replace(s$NAME, "_[1-4]$",""))

# load annotations
sub_annot_loc <- here("seqremoved_new")
bambu_sub_annot <- read.delim(here(sub_annot_loc, "/bambu/bambu_1.0.3_extended_annotations.gtf"), header = FALSE)
flair_sub_annot <- read.delim(here(sub_annot_loc,"/FLAIR/flair.collapse.isoforms.gtf"), header = FALSE)
flames_sub_annot <- read.delim(here(sub_annot_loc,"/FLAMES/isoform_annotated.filtered.gff3"), header = FALSE)
sqanti_sub_annot <- read.delim(here(sub_annot_loc, "/cupcake/merged_ont_sub.collapsed_corrected.gtf"), header = FALSE)
sqanti_class_sub <- read.delim(here(sub_annot_loc,"/cupcake/merged_ont_sub.collapsed_classification.txt"))
stringtie_sub_annot <- read.delim(here(sub_annot_loc,"/stringtie/stringtie_merged_sub_guided.gtf"), header = FALSE)
talon_sub_annot <- read.delim(here(sub_annot_loc, "/TALON/ONT_sub_filtered_annot_talon.gtf"), header = FALSE)

# load counts
sub_counts_loc <- here("seqremoved_new")
sub_bambu_counts <- read.delim(paste0(sub_counts_loc,"/bambu/bambu_1.0.3_counts_transcript.txt")) %>% data.frame(., total_A=rowSums(.[3:5]), total_B=rowSums(.[6:8]))
sub_flair_counts <- read.delim(paste0(sub_counts_loc,"/FLAIR/counts_matrix.tsv")) %>% data.frame(., total_A=rowSums(.[2:4]), total_B=rowSums(.[5:7]))
sub_flames_counts <- read_csv(paste0(sub_counts_loc,"/FLAMES/transcript_count.csv.gz")) %>% data.frame(., total_A=rowSums(.[3:5]), total_B=rowSums(.[6:8]))
sub_sqanti_counts <- read.delim(paste0(sub_counts_loc,"/ont_sub_cupcake/transcript_quant/quant.sf"))
sub_stringtie_counts <- read.delim(paste0(sub_counts_loc,"/ont_sub_stringtie/transcript_quant/quant.sf"))
sub_talon_counts <- read.delim(paste0(sub_counts_loc,"/TALON/ONT_sub_talon_abundance_filtered.tsv")) %>% data.frame(., total_A=rowSums(.[12:14]), total_B=rowSums(.[15:17]))

# isolate sequins-related transcript annotations (chrIS)
s_bambu_sub <- bambu_sub_annot[bambu_sub_annot$V1=="chrIS" & bambu_sub_annot$V3=="transcript",c(1:5,9)]
s_flair_sub <- flair_sub_annot[flair_sub_annot$V1=="chrIS" & flair_sub_annot$V3=="transcript",c(1:5,9)]
s_flames_sub <- flames_sub_annot[flames_sub_annot$V1=="chrIS" & flames_sub_annot$V3=="transcript",c(1:5,9)]
s_sqanti_sub <- sqanti_sub_annot[sqanti_sub_annot$V1=="chrIS" & sqanti_sub_annot$V3=="transcript",] %>% separate(., V9, into=c("isoform", "gene"), sep = ";")
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

summary_sub <- data.frame(tool = c("bambu", "Cupcake", "FLAIR", "FLAMES","StringTie2", "TALON"),
                          seq_related_annot = c(nrow(s_bambu_sub), nrow(s_sqanti_sub),nrow(s_flair_sub), nrow(s_flames_sub),nrow(s_stringtie_sub), nrow(s_talon_sub)),
                          known_seq_120 = c(nrow(s_known_bambu_sub), length(s_known_sqanti_sub),nrow(s_known_flair_sub), nrow(s_known_flames_sub), nrow(s_known_stringtie_sub), nrow(s_known_talon_sub)),
                          recovered_seq_40 = c(26,32,11,13,32,13))

summary_sub$kept_fd <- summary_sub$seq_related_annot - summary_sub$known_seq_120
summary_sub$removed_fd <- c(nrow(s_bambu_sub) - nrow(s_bambu_sub[str_detect(s_bambu_sub$V9,paste(removed_genes,collapse="|")),]),
                            nrow(s_sqanti_sub) - nrow(s_sqanti_sub[str_detect(s_sqanti_sub$associated_gene, paste(removed_genes,collapse="|")),]),
                            nrow(s_flair_sub) - nrow(s_flair_sub[str_detect(s_flair_sub$V9,paste(removed_genes,collapse="|")),]),
                            nrow(s_flames_sub) - nrow(s_flames_sub[str_detect(s_flames_sub$V9,paste(removed_genes,collapse="|")),]),
                            nrow(s_stringtie_sub) - nrow(s_stringtie_sub[str_detect(s_stringtie_sub$V9,paste(removed_genes,collapse="|")),]),
                            nrow(s_talon_sub) - nrow(s_talon_sub[str_detect(s_talon_sub$`1`,paste(removed_genes,collapse="|")),]))

summary_sub$precision_kept <- summary_sub$known_seq_120 / (summary_sub$known_seq_120 + summary_sub$kept_fd)
summary_sub$precision_removed <- summary_sub$recovered_seq_40 / (summary_sub$recovered_seq_40 + summary_sub$removed_fd)
summary_sub$recall_kept <- summary_sub$known_seq_120 / 120
summary_sub$recall_removed <- summary_sub$recovered_seq_40 / 40

data <- pivot_longer(summary_sub, c("precision_kept", "precision_removed", "recall_kept", "recall_removed"), names_to = "rate",values_to="value")
seq_sub <- ggbarplot(data, "tool", "value",
                     fill = "rate", color = "rate",
                     width =0.9,position = position_dodge(0.9),
                     palette = brewer.pal(5, "Set1")[c(2,3,4,5)],
                     xlab = "Method", ylab = "Rate",
                     facet.by="rate") +
  theme(legend.position="none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=15),
        axis.text.y = element_text(size=13), axis.title = element_text(size=15), strip.text.x = element_text(
          size = 13))


# pdf(here("sequins_recovered_rates_grid.pdf"),width=9, height=8)
seq_sub
# dev.off()



