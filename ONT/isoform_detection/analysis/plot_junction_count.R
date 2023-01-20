# Using short-read junction count from SQANTI3 to validate long-read detected isoforms
# Author: Xueyi Dong

library(ggplot2)
library(ggridges)
library(ggsci)
library(parallel)
library(plyr)
library(dplyr)

# prepare data----
# load("ont_counts_filtered.RData")
# junc_files <- list.files("./out", "junctions.txt", recursive = TRUE)
# junc <- lapply(junc_files,
#                function(x){
#                  junc_count = read.delim(file.path("./out", x))
#                  # method = strsplit2(x, "/")[1]
#                  # filt = junc_count$isoform %in% eval(parse(text = 
#                  #                                             paste0("ont_", method, "_counts")))[,1]
#                  return(junc_count)
#                })
# names(junc) <- limma::strsplit2(junc_files, "/")[,1]
# saveRDS(junc, "junc_unfiltered.RDS")
# lapply(junc, function(x){
#   length(unique(x$isoform))
# })
# # $bambu
# # [1] 261917
# # 
# # $flair
# # [1] 224055
# # 
# # $flames
# # [1] 75843
# # 
# # $sqanti
# # [1] 1292439
# # 
# # $stringtie2
# # [1] 204605
# # 
# # $talon
# # [1] 200800
# 
# # filter out lowly expressed transcripts
# table(junc$bambu$isoform %in% ont_bambu_counts$TXNAME)
# # FALSE   TRUE 
# # 631796 643252
# junc$bambu <- junc$bambu[junc$bambu$isoform %in% ont_bambu_counts$TXNAME, ]
# 
# flair_ids <- limma::strsplit2(ont_flair_counts$ids, "_")[ ,1]
# table(junc$flair$isoform %in% flair_ids)
# # FALSE    TRUE 
# # 736728 1016318
# junc$flair <- junc$flair[junc$flair$isoform %in% flair_ids, ]
# 
# flames_ids <- paste0("transcript:", ont_flames_counts$transcript_id)
# table(junc$flames$isoform %in% flames_ids)
# # TRUE 
# # 492562 
# junc$flames <- junc$flames[junc$flames$isoform %in% flames_ids, ]
# 
# table(junc$sqanti$isoform %in% ont_sqanti_counts$Name)
# # FALSE    TRUE 
# # 5550600 1037009 
# junc$sqanti <- junc$sqanti[junc$sqanti$isoform %in% ont_sqanti_counts$Name, ]
# 
# table(junc$stringtie2$isoform %in% ont_stringtie_counts$Name)
# # FALSE   TRUE 
# # 772705 428661
# junc$stringtie2 <- junc$stringtie2[junc$stringtie2$isoform %in% ont_stringtie_counts$Name, ]
# 
# table(junc$talon$isoform %in% ont_talon_counts$annot_transcript_id)
# # FALSE   TRUE 
# # 934563 123270 
# junc$talon <- junc$talon[junc$talon$isoform %in% ont_talon_counts$annot_transcript_id, ]
# saveRDS(junc, "junc.RDS")

# plot junction count distribution ---- 
junc <- readRDS(here("code","analysis","junc.RDS"))
lapply(junc, function(x){
  length(unique(x$isoform))
})
# deal with redundant junction information
junc_no_dup <- lapply(junc, function(x){
  x$junction = paste(x$chrom, x$genomic_start_coord, x$genomic_end_coord, sep = "_")
  x <- x[!duplicated(x$junction), ]
})
junc_df <- Reduce(rbind, junc_no_dup)
junc_df$Tool <- rep(c("bambu", "FLAIR", "FLAMES", "Cupcake", "StringTie2", "TALON"), 
                    sapply(junc_no_dup, nrow, simplify = TRUE))

# plot percentage of novel junctions supported by Illumina
# defind >10 counts as supported
junc_df$supported_by_Illumina <- junc_df$total_coverage_unique > 10
pdf("NovelJuncIlluminaSupport.pdf", height = 5, width = 8)
junc_df %>% 
  filter(junction_category == "novel") %>% 
  count(Tool = factor(Tool), supported_by_Illumina = factor(supported_by_Illumina)) %>% 
  plyr::ddply(., .(Tool), transform, pct = n / sum(n)) %>% 
  # mutate(pct = prop.table(n)) %>% 
  ggplot(aes(x = Tool, y = n, fill = supported_by_Illumina, label = scales::percent(pct))) +
  geom_bar(stat = "identity") +
  geom_text(position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("#A8D1D1", "#FD8A8A")) +
  labs(y = "Novel junctions", fill = "Supported by Illumina\n(counts > 10)") +
  theme_bw()+
  theme(text = element_text(size = 15))
dev.off()


ggplot(junc_df[junc_df$junction_category == "novel", ], aes(x = Tool, fill = supported_by_Illumina)) +
  geom_bar() +
  theme_bw() +
  theme(text = element_text(size = 15))


#pdf("plots/juncCovDistr.pdf", height = 4, width = 8)
junc_cov <- ggplot(junc_df, aes(x=(total_coverage_unique + 0.5), y=Tool, fill = Tool)) +
  geom_density_ridges(size = 0.3, scale = 3) +
  guides(fill=guide_legend(reverse = TRUE)) +
  scale_x_continuous(trans = "log10") +
  scale_y_discrete(expand=c(0,0))+
  facet_grid(rows = vars(junction_category)) +
  scale_fill_jama(alpha = 0.3) +
  labs(x="Splice junction Illumina read count") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_blank(),
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.line.x = element_line(size=0.1),axis.ticks.length=unit(0.25,"cm"), axis.text.x = element_text(size=12),
        text = element_text(size=15), legend.text = element_text(size=15)) 
#dev.off()
#pdf("plots/juncCovDistr2.pdf", height = 4, width = 4)
junc_cov_2 <- ggplot(junc_df, aes(x=(total_coverage_unique + 0.5), y=Tool, fill = Tool)) +
  geom_density_ridges(size = 0.3, scale = 3) +
  guides(fill=guide_legend(reverse = TRUE)) +
  scale_x_continuous(trans = "log10") +
  scale_y_discrete(expand=c(0,0))+
  scale_fill_jama(alpha = 0.3) +
  labs(x="Junction Illumina read count") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_blank(),
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.line.x = element_line(size=0.1),axis.ticks.length=unit(0.25,"cm"), axis.text.x = element_text(size=12),
        text = element_text(size=15), legend.text = element_text(size=15)) 
#dev.off()
# to do: add observation number


# with duplicated junctions
junc_df2 <- Reduce(rbind, junc)
junc_df2$Tool <- rep(c("bambu", "FLAIR", "FLAMES", "SQANTI3", "StringTie2", "TALON"), 
                    sapply(junc, nrow, simplify = TRUE))

# stats per isoform----
junc_df2$isoform_tool <- paste(junc_df2$isoform, junc_df2$Tool, sep = "_")
junc_df2$isDup <- duplicated(junc_df2$isoform_tool)
iso_stat <- data.frame(
  isoform = junc_df2$isoform[!junc_df2$isDup],
  tool = junc_df2$Tool[!junc_df2$isDup],
  isoform_tool = junc_df2$isoform_tool[!junc_df2$isDup]
)
# calculate number of junction per isoform
idx = which(!junc_df2$isDup)
idx2 = c(idx[-1], nrow(junc_df2)+1)
iso_stat$junc_count = idx2 - idx
junc_df2$supported = junc_df2$total_coverage_unique >= 10
iso_w_unsupported_junc = junc_df2$isoform_tool[!junc_df2$supported]
iso_stat$all_supported = !(iso_stat$isoform_tool %in% iso_w_unsupported_junc)
# calculate unsupported number for each junc
idx = which(!duplicated(iso_w_unsupported_junc))
idx2 = c(idx[-1], length(iso_w_unsupported_junc) + 1)
unsupport_info <- data.frame(
  isoform_tool = iso_w_unsupported_junc[!duplicated(iso_w_unsupported_junc)],
  number = idx2 - idx
)
iso_stat$junc_unsupport_count = 0
iso_stat$junc_unsupport_count[match(unsupport_info$isoform_tool, iso_stat$isoform_tool)] <- unsupport_info$number
iso_stat$junc_support_frac = 1 - (iso_stat$junc_unsupport_count / iso_stat$junc_count)

stat_box_data <- function(y, upper_limit = 1.15) {
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = length(y)
    )
  )
}
#pdf("plots/juncSupProp.pdf", height = 4, width = 6)
junc_sup_prop <- ggplot(iso_stat, aes(x=tool, y=junc_support_frac, fill=tool)) +
  geom_violin(alpha = 0.4) +
  geom_boxplot(width = 0.2, outlier.colour = NA, alpha = 0) +
  scale_fill_jama(alpha = 0.3) +
  theme_bw() + 
  labs(x = "Tool", y = "Proportion of supported junction") +
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text", 
    hjust = 0.5,
    vjust = 0.9
  )  +
  theme(text = element_text(size=18), legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))
#dev.off()

#pdf("plots/allJuncSup.pdf", height = 4, width = 6)
all_junc_sup <- ggplot(iso_stat, aes(x=tool, fill=all_supported)) +
  geom_bar() +
  labs(x = "Tool", y = "Number of transcripts", fill = "All junctions\nsupported by\nshort reads") +
  theme_bw() +
  theme(text = element_text(size=18), 
        axis.text.x = element_text(angle = 30, hjust = 1))
#dev.off()

# Look into isoform classification ----
class_files <- list.files("/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/ONT/isoform_detection/methods/junction_validation/out", "classification.txt", recursive = TRUE, full.names=TRUE)
class <- lapply(class_files,
               function(x){
                 class = read.delim(x)
                 return(class)
               })
names(class) <- limma::strsplit2(class_files, "/")[,1]

# remove FSM and remove redundant junction ----
junc_df3 = data.frame()
for(i in 1:6){
  tmp <- junc[[i]][class[[i]]$structural_category[match(junc[[i]]$isoform, class[[i]]$isoform)] != "full-splice_match", ]
  tmp$Tool = unique(junc_df2$Tool)[i]
  tmp$junction = paste(tmp$chrom, tmp$genomic_start_coord, tmp$genomic_end_coord, sep = "_")
  tmp <- tmp[!duplicated(tmp$junction), ]
  junc_df3 = rbind(junc_df3, tmp)
}
#pdf("plots/juncCovDistrNoFSM.pdf", height = 4, width = 8)
no_fsm <- ggplot(junc_df3, aes(x=(total_coverage_unique + 0.5), y=Tool, fill = Tool)) +
  geom_density_ridges(size = 0.3, scale = 3) +
  guides(fill=guide_legend(reverse = TRUE)) +
  scale_x_continuous(trans = "log10") +
  scale_y_discrete(expand=c(0,0))+
  facet_grid(rows = vars(junction_category)) +
  scale_fill_jama(alpha = 0.3) +
  labs(x="Splice junction Illumina read count") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_blank(),
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.line.x = element_line(size=0.1),axis.ticks.length=unit(0.25,"cm"), axis.text.x = element_text(size=12),
        text = element_text(size=15), legend.text = element_text(size=15)) 
#dev.off()
#pdf("plots/juncCovDistrNoFSM2.pdf", height = 4, width = 4)
no_fsm_2 <- ggplot(junc_df3, aes(x=(total_coverage_unique + 0.5), y=Tool, fill = Tool)) +
  geom_density_ridges(size = 0.3, scale = 3) +
  guides(fill=guide_legend(reverse = TRUE)) +
  scale_x_continuous(trans = "log10") +
  scale_y_discrete(expand=c(0,0))+
  scale_fill_jama(alpha = 0.3) +
  labs(x="Junction Illumina read count") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_blank(),
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.line.x = element_line(size=0.1),axis.ticks.length=unit(0.25,"cm"), axis.text.x = element_text(size=15),
        text = element_text(size=15), legend.text = element_text(size=15)) 
#dev.off()