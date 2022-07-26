cd /stornext/General/data/user_managed/grpu_mritchie_1/Mixture/MainStudy/in_silico_mix
for SAMPLE in R{1..3}
do
	cat subsample_fq/$SAMPLE-000.025.fastq.gz subsample_fq/$SAMPLE-100.075.fastq.gz > $SAMPLE-075sub.fastq.gz
	cat subsample_fq/$SAMPLE-000.050.fastq.gz subsample_fq/$SAMPLE-100.050.fastq.gz > $SAMPLE-050sub.fastq.gz
	cat subsample_fq/$SAMPLE-000.075.fastq.gz subsample_fq/$SAMPLE-100.025.fastq.gz > $SAMPLE-025sub.fastq.gz
done