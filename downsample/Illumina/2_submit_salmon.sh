for SAMPLE in H1975-1_S9  H1975-2_S10  H1975-3_S11  HCC827-1_S12 HCC827-2_S13_merged HCC827-5-repeat_S14 mix_025075_{1..3} mix_050050_{1..3} mix_075025_{1..3}
do
	sbatch 2_salmon.sh $SAMPLE 1M
	sbatch 2_salmon.sh $SAMPLE 3M
	sbatch 2_salmon.sh $SAMPLE 5M
	sbatch 2_salmon.sh $SAMPLE 10M
	sbatch 2_salmon.sh $SAMPLE 15M
	sbatch 2_salmon.sh $SAMPLE 20M
done