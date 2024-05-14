#!/bin/bash -l
#SBATCH --partition=epyc
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --mem=100g
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name="02-rnaseq-diffexpr"
#SBATCH --output=log/%x_%A_%a.log
##################################################################

## RNA-seq Workflow - Alignment #############
##
## This workflow will carry out the following steps:
## 1 - RNA-seq differential abundance
## 2 - post-analysis QC
## 3 - Generate DEG lists
## 4 - Summary plots

## Get configuration file
if [ -f code/config.txt ]; then
  source code/config.txt
else
    echo "You're missing the config.txt file"
    exit
fi

# get starting time point
echo "Run start"
timestamp

## Job variables
N=${SLURM_ARRAY_TASK_ID}
NUM_CORES=8

## Generate output directories
[ ! -d "$DESEQ2_DIR" ] && \
mkdir -p "${DESEQ2_DIR}/genelists" && \
mkdir -p "${DESEQ2_DIR}/plots"

## Running workflow
echo "## Starting Differential abundance analysis ##########"

Rscript --vanilla code/DEG_analysis.R ${DESEQ2_FCT}
multiqc_summary

# end time stamp
echo "Run complete"
timestamp