library(data.table)

Files <-
  "/stornext/Projects/promethion/promethion_access/lab_ritchie/transcr_bench/long_term/transcr_bench/guppy_4_rebasecall_Oct20/sequencing_summary.txt.gz"
Table <- fread(Files, header = TRUE, sep = "\t")
Read_Id <- as.character(Table$read_id)
Read_length <- as.numeric(Table$sequence_length_template)
Qscore <- as.numeric(Table$mean_qscore_template)
Pass_filtering <- Table$pass_filtering
Barcode <- as.character(Table$barcode_arrangement)
Table <- cbind(Read_Id, Read_length, Qscore, Pass_filtering, Barcode)
saveRDS(Table, "summaryInfo.RDS")