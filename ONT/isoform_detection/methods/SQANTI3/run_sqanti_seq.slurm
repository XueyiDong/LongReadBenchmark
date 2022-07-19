#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --output=sqanti_ont.out
#SBATCH --mail-user=du.m@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

# Author: Mei Du
# Using SQANTI3 version 1.6.0 (https://github.com/ConesaLab/SQANTI3)
# Run on the full dataset with the modified GTF

module load anaconda3

COLLAPSED=merged_ont.collapsed.gff
ANNOT=/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/references/genome_rnasequin_decoychr_2.4_edited_seqremoved.gtf
GENOME=/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/references/genome_rnasequin_decoychr_2.4.fa

source activate /stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/SQANTI3/SQANTI3_env
export PYTHONPATH=$PYTHONPATH:/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/SQANTI3/cDNA_Cupcake/sequence
export PYTHONPATH=$PYTHONPATH:/stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/SQANTI3/cDNA_Cupcake/

python /stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/SQANTI3/sqanti3_qc.py \
--gtf $COLLAPSED \
$ANNOT \
$GENOME > sqanti_ont.out
