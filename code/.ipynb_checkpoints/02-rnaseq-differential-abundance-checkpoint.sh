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
## 1 - pre-analysis QC
## 2 - RNA-seq differential abundance
## 3 - post-analysis QC
## 4 - Generate DEG lists

## Get configuration file
if [ -f code/config.txt ]; then
  source code/config.txt
else
    echo "You're missing the config.txt file"
    exit
fi

## Load modules and environments
module load ${MODULE_LIST}

N=${SLURM_ARRAY_TASK_ID}
NUM_CORES=8

## Generate output directories

[ ! -d "$DESEQ2_DIR" ] && \
mkdir -p "${$DESEQ2_DIR}/genelists" && \
mkdir -p "${$DESEQ2_DIR}/plots"

## Running workflow

echo "## Starting Differential abundance analysis ##########"

sed -n ${N}p $METADATA | while IFS="," read 

do

  # Rscript --vanilla DEG_pre-processing.R
  Rscript --vanilla DEG_analysis.R
  # Rscript --vanilla DEG_post-processing.R

done

