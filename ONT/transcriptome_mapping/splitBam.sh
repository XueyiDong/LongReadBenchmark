#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL
module load samtools
BC=$1
samtools view -H barcode$BC.primary.bam > header$BC
samtools view barcode$BC.primary.bam | split - "./split/barcode$BC" -l 8000000 --filter='cat header$BC - | samtools view -b - > $FILE.bam'