FQ_PURE=/stornext/Projects/promethion/promethion_access/lab_ritchie/transcr_bench/long_term/transcr_bench/guppy_4_rebasecall_Oct20/barcode01-06
OUT=/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/downsample/ONT/fastq
for SAMPLE in barcode{01..06}
do
	sbatch 1_downsamplefq.sh $SAMPLE 1000000 $OUT/$SAMPLE.1M $FQ_PURE
	sbatch 1_downsamplefq.sh $SAMPLE 3000000 $OUT/$SAMPLE.3M $FQ_PURE
	sbatch 1_downsamplefq.sh $SAMPLE 5000000 $OUT/$SAMPLE.5M $FQ_PURE
	sbatch 1_downsamplefq.sh $SAMPLE 10000000 $OUT/$SAMPLE.10M $FQ_PURE
	sbatch 1_downsamplefq.sh $SAMPLE 15000000 $OUT/$SAMPLE.15M $FQ_PURE
	sbatch 1_downsamplefq.sh $SAMPLE 20000000 $OUT/$SAMPLE.20M $FQ_PURE
done

FQ_MIX=/stornext/General/data/user_managed/grpu_mritchie_0/XueyiDong/long_read_benchmark/ONT/mix_fq/
for SAMPLE in mix_025075_{1..3} mix_050050_{1..3} mix_075025_{1..3}
do
	sbatch 1_downsamplefq.sh $SAMPLE 1000000 $OUT/$SAMPLE.1M $FQ_MIX
	sbatch 1_downsamplefq.sh $SAMPLE 5000000 $OUT/$SAMPLE.5M $FQ_MIX
	sbatch 1_downsamplefq.sh $SAMPLE 10000000 $OUT/$SAMPLE.10M $FQ_MIX
	sbatch 1_downsamplefq.sh $SAMPLE 15000000 $OUT/$SAMPLE.15M $FQ_MIX
	sbatch 1_downsamplefq.sh $SAMPLE 20000000 $OUT/$SAMPLE.20M $FQ_MIX
done