cd /stornext/General/data/user_managed/grpu_mritchie_1/Mixture/MainStudy/in_silico_mix
for SAMPLE in R{1..3}
do
	cat subsample_fq/$SAMPLE-000.025.fastq.gz subsample_fq/$SAMPLE-100.075.fastq.gz > $SAMPLE-075sub.fastq.gz
	cat subsample_fq/$SAMPLE-000.050.fastq.gz subsample_fq/$SAMPLE-100.050.fastq.gz > $SAMPLE-050sub.fastq.gz
	cat subsample_fq/$SAMPLE-000.075.fastq.gz subsample_fq/$SAMPLE-100.025.fastq.gz > $SAMPLE-025sub.fastq.gz
done

# cat R1-000-sub25.fastq R1-100-sub75.fastq > R1-075sub.fastq
# cat R1-000-sub50.fastq R1-100-sub50.fastq > R1-050sub.fastq
# cat R1-000-sub75.fastq R1-100-sub25.fastq > R1-025sub.fastq
# cat R2-000-sub25.fastq R2-100-sub75.fastq > R2-075sub.fastq
# cat R2-000-sub50.fastq R2-100-sub50.fastq > R2-050sub.fastq
# cat R2-000-sub75.fastq R2-100-sub25.fastq > R2-025sub.fastq
# cat R3-000-sub25.fastq R3-100-sub75.fastq > R3-075sub.fastq
# cat R3-000-sub50.fastq R3-100-sub50.fastq > R3-050sub.fastq
# cat R3-000-sub75.fastq R3-100-sub25.fastq > R3-025sub.fastq