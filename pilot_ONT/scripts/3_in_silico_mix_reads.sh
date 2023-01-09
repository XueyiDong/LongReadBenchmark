cd /vast/scratch/users/dong.x/long_read_benchmark/pilot_ONT/
# for SAMPLE in R{1..3}
# do
# 	cat subsample_fq/$SAMPLE-000.025.fastq.gz subsample_fq/$SAMPLE-100.075.fastq.gz > $SAMPLE-075sub.fastq.gz
# 	cat subsample_fq/$SAMPLE-000.050.fastq.gz subsample_fq/$SAMPLE-100.050.fastq.gz > $SAMPLE-050sub.fastq.gz
# 	cat subsample_fq/$SAMPLE-000.075.fastq.gz subsample_fq/$SAMPLE-100.025.fastq.gz > $SAMPLE-025sub.fastq.gz
# done

mkdir -p fq_mix
cat fq_subsample/barcode01.025.fq.gz fq_subsample/barcode05.075.fq.gz > fq_mix/rep1-075.fq.gz
cat fq_subsample/barcode01.050.fq.gz fq_subsample/barcode05.050.fq.gz > fq_mix/rep1-050.fq.gz
cat fq_subsample/barcode01.075.fq.gz fq_subsample/barcode05.025.fq.gz > fq_mix/rep1-025.fq.gz
cat fq_subsample/barcode06.025.fq.gz fq_subsample/barcode10.075.fq.gz > fq_mix/rep2-075.fq.gz
cat fq_subsample/barcode06.050.fq.gz fq_subsample/barcode10.050.fq.gz > fq_mix/rep2-050.fq.gz
cat fq_subsample/barcode06.075.fq.gz fq_subsample/barcode10.025.fq.gz > fq_mix/rep2-025.fq.gz
