WD=/vast/scratch/users/dong.x/long_read_benchmark/pilot_ONT

for SAMPLE in rep1-025 rep1-050 rep1-075 rep2-025 rep2-050 rep2-075
do
	mkdir -p $WD/salmon/$SAMPLE
	sbatch 5_salmon_quant.sh $WD/salmon/$SAMPLE $WD/bam_tx/$SAMPLE.sorted
done