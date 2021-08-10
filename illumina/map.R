library(Rsubread)
datadir="/stornext/General/data/user_managed/grpu_mritchie_1/long_RNA_benchmark/Illumina_April20/data/raw/"
index="/wehisan/home/allstaff/d/dong.x/annotation/HumanSequins/GrCh38_sequins"
outdir="/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/illumina/bam/"
# 
# subread-align -t 0 -T 8 -i $index -r $datadir/H1975-1_S9_R1.fastq.gz -R $datadir/H1975-1_S9_R2.fastq.gz -o $outdir/H1975-1_S9.bam
# subread-align -t 0 -T 8 -i $index -r $datadir/H1975-2_S10_R1.fastq.gz -R $datadir/H1975-2_S10_R2.fastq.gz -o $outdir/H1975-2_S10.bam
# subread-align -t 0 -T 8 -i $index -r $datadir/H1975-3_S11_R1.fastq.gz -R $datadir/H1975-3_S11_R2.fastq.gz -o $outdir/H1975-3_S11.bam
# subread-align -t 0 -T 8 -i $index -r $datadir/HCC827-1_S12_R1.fastq.gz -R $datadir/HCC827-1_S12_R2.fastq.gz -o $outdir/HCC827-1_S12.bam
# subread-align -t 0 -T 8 -i $index -r $datadir/HCC827-2_S13_R1.fastq.gz -R $datadir/HCC827-2_S13_R2.fastq.gz -o $outdir/HCC827-2_S13.bam
# subread-align -t 0 -T 8 -i $index -r $datadir/HCC827-5-repeat_S14_R1.fastq.gz -R $datadir/HCC827-5-repeat_S14_R2.fastq.gz -o $outdir/HCC827-5-repeat_S14.bam

alignIllumina <- function(prefix){
  align(index = index, 
        readfile1 = paste0(datadir, prefix, "_R1.fastq.gz"),
        readfile2 = paste0(datadir, prefix, "_R2.fastq.gz"),
        output_file = paste0(outdir, prefix, ".bam"),
        nthreads = 8)
}

alignIllumina("H1975-1_S9")
alignIllumina("H1975-2_S10")
alignIllumina("H1975-3_S11")
alignIllumina("HCC827-1_S12")
alignIllumina("HCC827-2_S13")
alignIllumina("HCC827-2_S13_topup")
alignIllumina("HCC827-5-repeat_S14")

files = file.path(outdir, list.files(outdir, ".bam$"))
files <- files[c(1, 3, 5, 7, 10, 12, 14)]
fc <- featureCounts(files = files,
                    annot.ext = "/wehisan/home/allstaff/d/dong.x/annotation/HumanSequins/gencode.v33.sequins.gtf",
                    isGTFAnnotationFile = TRUE,
                    isPairedEnd = TRUE,
                    nthreads = 8
                    )
saveRDS(fc, file="counts.RDS")
