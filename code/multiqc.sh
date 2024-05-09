#!/bin/bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=12g
#SBATCH --time=0-2:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name="multiqc"
#SBATCH --output=log/%x_%j.log

##################################################################
## Shell script to generate MultiQC summary

## Load environment
source activate DEG-analysis

## Setting variables
OUT_DIR=analysis/multiqc

# Create analysis folder (if it doesn't exist)
[ ! -d "$OUT_DIR" ] && mkdir -p "$OUT_DIR"

# Run multiqc
multiqc --force \
    --dirs \
    --outdir ${OUT_DIR} \
    --filename multiqc \
    analysis