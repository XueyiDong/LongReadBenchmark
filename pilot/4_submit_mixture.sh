for SAMPLE in R{1..3}
# make 40M reads in total
do
	sbatch 4_map_mixture_samples.sh $SAMPLE-025sub
	sbatch 4_map_mixture_samples.sh $SAMPLE-050sub
	sbatch 4_map_mixture_samples.sh $SAMPLE-075sub
done