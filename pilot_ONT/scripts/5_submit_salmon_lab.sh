WD=/vast/scratch/users/dong.x/long_read_benchmark/pilot_ONT

for SAMPLE in barcode{01..12} unclassified
do
	mkdir -p $WD/salmon/$SAMPLE\_pass
	sbatch 5_salmon_quant.sh $WD/salmon/$SAMPLE\_pass $WD/bam_tx/$SAMPLE\_sup_q10_pass.sorted
done