#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --job-name=SQANTI3
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

# Author: Mei Du, Xueyi Dong
# Using SQANTI3 version 1.6.0 (https://github.com/ConesaLab/SQANTI3)

module load anaconda3

OUT_DIR="/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/ONT/isoform_detection/methods/junction_validation"

cd ${OUT_DIR}
source activate SQANTI3.env
mkdir -p out/$2
cd out/$2
export PYTHONPATH=$PYTHONPATH:~/Programs/cDNA_Cupcake/sequence/
export PYTHONPATH=$PYTHONPATH:~/Programs/cDNA_Cupcake/
echo $2
python ~/Programs/SQANTI3-5.0/sqanti3_qc.py \
/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/plots/readcount_cor/ONT_annotations/$1 \
/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/references/genome_rnasequin_decoychr_2.4_edited.gtf \
/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/references/genome_rnasequin_decoychr_2.4.fa \
--cpus 16 --skipORF --short_reads ../../shortRead.fofn > sqanti_$2.out
# --SR_bam ../../shortBam.fofn \
