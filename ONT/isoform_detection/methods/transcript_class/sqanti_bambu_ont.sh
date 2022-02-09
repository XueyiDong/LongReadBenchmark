# Author: Mei Du
# Using SQANTI3 version 1.6.0 (https://github.com/ConesaLab/SQANTI3)

module load anaconda3

OUT_DIR="/stornext/General/data/user_managed/grpu_mritchie_0/Mei/long_read_benchmark/SQANTI3/bambu"
GTF="/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/bambu/outputs/ONT_extended_annotations.gtf"

cd ${OUT_DIR}
source activate /stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/SQANTI3/SQANTI3_env
export PYTHONPATH=$PYTHONPATH:/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/SQANTI3/cDNA_Cupcake/sequence

python /stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/SQANTI3/sqanti3_qc.py \
--gtf ${GTF} \
/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/references/genome_rnasequin_decoychr_2.4_edited.gtf \
/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/references/genome_rnasequin_decoychr_2.4.fa > sqanti_bambu_ont.out

