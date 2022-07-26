for SAMPLE in H1975-1_S9  H1975-2_S10  H1975-3_S11  HCC827-1_S12 HCC827-2_S13_merged HCC827-5-repeat_S14
do
	sbatch downsamplefq.sh $SAMPLE\_R1 7500000 $SAMPLE\_R1.025
	sbatch downsamplefq.sh $SAMPLE\_R1 15000000 $SAMPLE\_R1.050
	sbatch downsamplefq.sh $SAMPLE\_R1 22500000 $SAMPLE\_R1.075
	sbatch downsamplefq.sh $SAMPLE\_R2 7500000 $SAMPLE\_R2.025
	sbatch downsamplefq.sh $SAMPLE\_R2 15000000 $SAMPLE\_R2.050
	sbatch downsamplefq.sh $SAMPLE\_R2 22500000 $SAMPLE\_R2.075
done