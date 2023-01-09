fq=/stornext/Projects/promethion/promethion_access/lab_ritchie/Xueyi_cellLine_mixture/long_term/fastq/H1_1-HC_5/20221107_1451_2D_PAK15503_0f50058f/guppy6.2.1_sup_prom
for SAMPLE in barcode{01..12} unclassified
do
	wc -l  $fq/$SAMPLE\_sup_q10_pass
done