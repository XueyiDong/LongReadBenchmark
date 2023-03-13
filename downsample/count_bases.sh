#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=count_bases
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --output=count_bases.out
#SBATCH --mail-type=FAIL

module load stornext
echo $1
snretrieve $1
zcat $1 | awk '{if (NR%4 == 2) print $0; }' | wc -cl > $2
snrmdiskcopy $1