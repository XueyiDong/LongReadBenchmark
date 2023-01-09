fq=/vast/scratch/users/dong.x/long_read_benchmark/pilot_ONT/fq_mix
for SAMPLE in rep1-025 rep1-050 rep1-075 rep2-025 rep2-050 rep2-075
do
	sbatch 4_map_ONT_tx.sh $fq $SAMPLE
done