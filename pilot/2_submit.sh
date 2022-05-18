for SAMPLE in R{1..3}-100 R{1..3}-000
# make 40M reads in total
do
	sbatch 2_downsamplefq.sh $SAMPLE 10000000 $SAMPLE.025
	sbatch 2_downsamplefq.sh $SAMPLE 20000000 $SAMPLE.050
	sbatch 2_downsamplefq.sh $SAMPLE 30000000 $SAMPLE.075
done