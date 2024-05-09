#!/bin/bash -l
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=16
#SBATCH --mem=100g
#SBATCH --time=0-8:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name="deseq2"
#SBATCH --output=log/%x_%j.log
##################################################################

####################
# Executing Rscript
# conda activate DEG-analysis
# Rscript --vanilla DESeq2.R
####################

# Load conda environment
source activate DEG-analysis

# Load variables
GENELISTS_DIR=analysis/deseq2/genelists
PLOTS_DIR=analysis/deseq2/plots

# Create analysis folder (if it doesn't exist)
[ ! -d "$GENELISTS_DIR" ] && mkdir -p "$GENELISTS_DIR"
[ ! -d "$PLOTS_DIR" ] && mkdir -p "$PLOTS_DIR"

# Run R script
Rscript --vanilla code/DESeq2.R