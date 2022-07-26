for SAMPLE in mix_025075_{1..3} mix_050050_{1..3} mix_075025_{1..3}
do
	sbatch salmon.sh $SAMPLE
done