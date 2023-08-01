WD=/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/downsample/ONT
for SAMPLE in `ls /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/downsample/ONT/bam | grep "barcode"`
do
	mkdir -p $WD/salmon/${SAMPLE%.bam*}
	sbatch 3_salmon_quant.sh $WD/salmon/${SAMPLE%.bam*} $WD/bam/${SAMPLE%.bam*}
done