#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --output=star_index.out
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

module load STAR/2.7.9a

cd /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/HumanSequins/

STAR --runThreadN 16 --runMode genomeGenerate --genomeFastaFiles GrCh38_sequins.fa \
--sjdbGTFfile gencode.v33.sequins.gtf --sjdbOverhang 80