#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --cpus-per-task=32
#SBATCH --output=log_fastqc.txt

for file in ./*.fastq; do
srun fastqc "$file"
done
