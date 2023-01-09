for SAMPLE in barcode01 barcode05 barcode06 barcode10
# make 3M reads in total each sample
do
	sbatch 2_downsamplefq.sh $SAMPLE 750000 $SAMPLE.025
	sbatch 2_downsamplefq.sh $SAMPLE 1500000 $SAMPLE.050
	sbatch 2_downsamplefq.sh $SAMPLE 2250000 $SAMPLE.075
done