OUT=/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/downsample/Illumina/fastq

FQ_PURE=/stornext/General/data/user_managed/grpu_mritchie_1/long_RNA_benchmark/Illumina_April20/data/raw/
for SAMPLE in H1975-1_S9  H1975-2_S10  H1975-3_S11  HCC827-1_S12 HCC827-5-repeat_S14
do
	sbatch 1_downsamplefq.sh $SAMPLE\_R1 1000000 $OUT/$SAMPLE\_R1.1M $FQ_PURE
	sbatch 1_downsamplefq.sh $SAMPLE\_R1 3000000 $OUT/$SAMPLE\_R1.3M $FQ_PURE
	sbatch 1_downsamplefq.sh $SAMPLE\_R1 5000000 $OUT/$SAMPLE\_R1.5M $FQ_PURE
	sbatch 1_downsamplefq.sh $SAMPLE\_R1 10000000 $OUT/$SAMPLE\_R1.10M $FQ_PURE
	sbatch 1_downsamplefq.sh $SAMPLE\_R1 15000000 $OUT/$SAMPLE\_R1.15M $FQ_PURE
	sbatch 1_downsamplefq.sh $SAMPLE\_R1 20000000 $OUT/$SAMPLE\_R1.20M $FQ_PURE
	sbatch 1_downsamplefq.sh $SAMPLE\_R2 1000000 $OUT/$SAMPLE\_R2.1M $FQ_PURE
	sbatch 1_downsamplefq.sh $SAMPLE\_R2 1000000 $OUT/$SAMPLE\_R2.1M $FQ_PURE
	sbatch 1_downsamplefq.sh $SAMPLE\_R2 3000000 $OUT/$SAMPLE\_R2.3M $FQ_PURE
	sbatch 1_downsamplefq.sh $SAMPLE\_R2 5000000 $OUT/$SAMPLE\_R2.5M $FQ_PURE
	sbatch 1_downsamplefq.sh $SAMPLE\_R2 10000000 $OUT/$SAMPLE\_R2.10M $FQ_PURE
	sbatch 1_downsamplefq.sh $SAMPLE\_R2 15000000 $OUT/$SAMPLE\_R2.15M $FQ_PURE
	sbatch 1_downsamplefq.sh $SAMPLE\_R2 20000000 $OUT/$SAMPLE\_R2.20M $FQ_PURE
done

FQ_S13=/stornext/General/data/user_managed/grpu_mritchie_0/XueyiDong/long_read_benchmark/illumina/
for SAMPLE in HCC827-2_S13_merged
do
	sbatch 1_downsamplefq.sh $SAMPLE\_R1 1000000 $OUT/$SAMPLE\_R1.1M $FQ_S13
	sbatch 1_downsamplefq.sh $SAMPLE\_R1 3000000 $OUT/$SAMPLE\_R1.3M $FQ_S13
	sbatch 1_downsamplefq.sh $SAMPLE\_R1 5000000 $OUT/$SAMPLE\_R1.5M $FQ_S13
	sbatch 1_downsamplefq.sh $SAMPLE\_R1 10000000 $OUT/$SAMPLE\_R1.10M $FQ_S13
	sbatch 1_downsamplefq.sh $SAMPLE\_R1 15000000 $OUT/$SAMPLE\_R1.15M $FQ_S13
	sbatch 1_downsamplefq.sh $SAMPLE\_R1 20000000 $OUT/$SAMPLE\_R1.20M $FQ_S13
	sbatch 1_downsamplefq.sh $SAMPLE\_R2 1000000 $OUT/$SAMPLE\_R2.1M $FQ_S13
	sbatch 1_downsamplefq.sh $SAMPLE\_R2 3000000 $OUT/$SAMPLE\_R2.3M $FQ_S13
	sbatch 1_downsamplefq.sh $SAMPLE\_R2 5000000 $OUT/$SAMPLE\_R2.5M $FQ_S13
	sbatch 1_downsamplefq.sh $SAMPLE\_R2 10000000 $OUT/$SAMPLE\_R2.10M $FQ_S13
	sbatch 1_downsamplefq.sh $SAMPLE\_R2 15000000 $OUT/$SAMPLE\_R2.15M $FQ_S13
	sbatch 1_downsamplefq.sh $SAMPLE\_R2 20000000 $OUT/$SAMPLE\_R2.20M $FQ_S13
done

FQ_MIX=/stornext/General/data/user_managed/grpu_mritchie_0/XueyiDong/long_read_benchmark/illumina/mix_fq
for SAMPLE in mix_025075_{1..3} mix_050050_{1..3} mix_075025_{1..3}
do
	sbatch 1_downsamplefq2.sh $SAMPLE\_R1 1000000 $OUT/$SAMPLE\_R1.1M $FQ_MIX
	sbatch 1_downsamplefq2.sh $SAMPLE\_R1 5000000 $OUT/$SAMPLE\_R1.5M $FQ_MIX
	sbatch 1_downsamplefq2.sh $SAMPLE\_R1 10000000 $OUT/$SAMPLE\_R1.10M $FQ_MIX
	sbatch 1_downsamplefq2.sh $SAMPLE\_R1 15000000 $OUT/$SAMPLE\_R1.15M $FQ_MIX
	sbatch 1_downsamplefq2.sh $SAMPLE\_R1 20000000 $OUT/$SAMPLE\_R1.20M $FQ_MIX
	sbatch 1_downsamplefq2.sh $SAMPLE\_R2 1000000 $OUT/$SAMPLE\_R2.1M $FQ_MIX
	sbatch 1_downsamplefq2.sh $SAMPLE\_R2 5000000 $OUT/$SAMPLE\_R2.5M $FQ_MIX
	sbatch 1_downsamplefq2.sh $SAMPLE\_R2 10000000 $OUT/$SAMPLE\_R2.10M $FQ_MIX
	sbatch 1_downsamplefq2.sh $SAMPLE\_R2 15000000 $OUT/$SAMPLE\_R2.15M $FQ_MIX
	sbatch 1_downsamplefq2.sh $SAMPLE\_R2 15000000 $OUT/$SAMPLE\_R2.20M $FQ_MIX
done