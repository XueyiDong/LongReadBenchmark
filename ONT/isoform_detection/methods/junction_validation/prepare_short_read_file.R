dir = "/stornext/General/data/user_managed/grpu_mritchie_1/long_RNA_benchmark/Illumina_April20/data/raw"
short_fq = list.files(dir, "fastq.gz")
R1 <- short_fq[grepl("R1", short_fq)]
R2 <- short_fq[grepl("R2", short_fq)]
short_fq = data.frame(
  R1 = file.path(dir, R1),
  R2 = file.path(dir, R2)
)
write.table(short_fq, "shortRead.fofn",
            row.names=FALSE,
            col.names=FALSE,
            sep = " ",
            quote = FALSE)

bam.dir = "/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/illumina/bam"
bam.sorted = list.files(bam.dir, "sorted.bam$")
bam = data.frame(
  bam = file.path(bam.dir, bam.sorted[c(-6, -7)])
)
write.table(bam, "shortBam.fofn",
            row.names=FALSE,
            col.names=FALSE,
            sep = " ",
            quote = FALSE)
