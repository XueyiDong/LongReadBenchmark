# Author: Mei Du
# Using TranscriptClean version 2.0.2 (https://github.com/mortazavilab/TranscriptClean)

# trimmed FASTA headers first after TranscriptClean error
# cut -d ' ' -f1 genome_rnasequin_decoychr_2.4.fa > genome_rnasequin_decoychr_2.4_trimmed.fa
cd /stornext/General/data/user_managed/grpu_mritchie_0/Mei/long_read_benchmark/TALON
for FILE in *.sam; do
touch "transcriptclean_ont_${FILE%'_sorted.sam'*}.slurm"
slurm_file="transcriptclean_ont_${FILE%'_sorted.sam'*}.slurm"

echo "#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --job-name=tc
#SBATCH --output=transcriptclean_ont_${FILE%'_sorted.sam'*}.out
#SBATCH --mail-user=du.m@wehi.edu.au
#SBATCH --mail-type=END,FAIL

module load anaconda3
source activate
conda activate TranscriptClean

python /stornext/General/data/user_managed/grpu_mritchie_0/Mei/long_read_benchmark/TALON/TranscriptClean-2.0.2/TranscriptClean.py \
--sam /stornext/General/data/user_managed/grpu_mritchie_0/Mei/long_read_benchmark/TALON/$FILE \
--genome /stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/references/genome_rnasequin_decoychr_2.4_trimmed.fa \
--outprefix ${FILE%'_sorted.sam'*} \
-t 16" > $slurm_file

done