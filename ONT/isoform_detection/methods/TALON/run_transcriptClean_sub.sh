# Author: Mei Du
# Using TranscriptClean version 2.0.2 (https://github.com/mortazavilab/TranscriptClean)

# trimmed FASTA headers first after TranscriptClean error
# cut -d ' ' -f1 genome_rnasequin_decoychr_2.4.fa > genome_rnasequin_decoychr_2.4_trimmed.fa
OUT=/vast/scratch/users/du.m/subsample/TALON
cd /stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/ONT_aligned/subsample_aligned
for FILE in barcode0[1-6]*.sam; do
touch $OUT/transcriptclean_ont_${FILE%'.sam'*}.slurm
slurm_file=$OUT/transcriptclean_ont_${FILE%'.sam'*}.slurm

echo "#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --job-name=tc
#SBATCH --output=transcriptclean_ont_${FILE%'.sam'*}.out
#SBATCH --mail-user=du.m@wehi.edu.au
#SBATCH --mail-type=END,FAIL

mkdir $OUT/${FILE%'.sam'*}
cd $OUT/${FILE%'.sam'*}

module load anaconda3
source activate
conda activate TranscriptClean

python /stornext/General/data/user_managed/grpu_mritchie_0/Mei/long_read_benchmark/TALON/TranscriptClean-2.0.2/TranscriptClean.py \
--sam /stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/ONT_aligned/subsample_aligned/$FILE \
--genome /stornext/General/data/user_managed/grpu_mritchie_1/Mei/long_read_benchmark/references/genome_rnasequin_decoychr_2.4_trimmed.fa \
--outprefix ONT_sub_TC \
--tmpDir $OUT/${FILE%'.sam'*} -t 8" > $slurm_file

sbatch $slurm_file

done