library(data.table)

Files <-
  "/stornext/Projects/promethion/promethion_access/lab_ritchie/Xueyi_cellLine_mixture/long_term/fastq/H1_1-HC_5/20221107_1451_2D_PAK15503_0f50058f/guppy6.2.1_sup_prom/sequencing_summary.txt.gz"
Table <- fread(Files, header = TRUE, sep = "\t")
Read_Id <- as.character(Table$read_id)
Read_length <- as.numeric(Table$sequence_length_template)
Qscore <- as.numeric(Table$mean_qscore_template)
Pass_filtering <- Table$passes_filtering
Barcode <- as.character(Table$barcode_arrangement)
Table <- cbind(Read_Id, Read_length, Qscore, Pass_filtering, Barcode)
saveRDS(Table, "summaryInfo.RDS")