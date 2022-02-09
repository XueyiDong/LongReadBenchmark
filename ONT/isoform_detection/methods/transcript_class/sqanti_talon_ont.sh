#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --output=sqanti_talon_ont.out
#SBATCH --mail-user=du.m@wehi.edu.au
#SBATCH --mail-type=END,FAIL

# Author: Mei Du
# Using SQANTI3 version 1.6.0 (https://github.com/ConesaLab/SQANTI3)

module load anaconda3

GTF="/stornext/General/data/user_managed/grpu_mritchie_0/Mei/long_read_benchmark/TALON/ONT/ONT_filtered_annot_talon_edited.gtf"

source activate /stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/SQANTI3/SQANTI3_env
export PYTHONPATH=$PYTHONPATH:/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/SQANTI3/cDNA_Cupcake/sequence

python /stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/SQANTI3/sqanti3_qc.py \
--gtf ${GTF} \
/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/references/genome_rnasequin_decoychr_2.4_edited.gtf \
/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/references/genome_rnasequin_decoychr_2.4.fa \
--skipORF

