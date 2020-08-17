#!/bin/bash
#SBATCH --job-name=salmon
#SBATCH --cpus-per-task=32
#SBATCH --output=log_salmon.txt

srun module load salmon/0.11.3
for r1_file in *_R1_*.fastq.gz; do
  r2_file="${r1_file/_R1_/_R2_}"
  out="${r1_file/_R1_/_both_}"
srun salmon quant -i /home/jrandelovic/data/Genomes/P_lividus/str_k25/ --gcBias -l IU -1 "$r1_file" -2 "$r2_file" -o mapped/quant_"$out"
done
