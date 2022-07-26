for SAMPLE in barcode{01..06}
do
	sbatch downsamplefq.sh $SAMPLE 7500000 $SAMPLE.025
	sbatch downsamplefq.sh $SAMPLE 15000000 $SAMPLE.050
	sbatch downsamplefq.sh $SAMPLE 22500000 $SAMPLE.075
done