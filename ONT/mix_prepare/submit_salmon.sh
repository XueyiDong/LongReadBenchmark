for SAMPLE in mix_025075_3 mix_050050_{1..3} mix_075025_{1..2}
do
	sbatch salmon_quant_mix.sh $SAMPLE
done