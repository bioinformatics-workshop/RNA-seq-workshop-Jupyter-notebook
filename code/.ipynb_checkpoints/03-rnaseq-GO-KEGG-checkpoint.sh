#!/bin/bash -l
#SBATCH --partition=epyc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100g
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name="03-rnaseq-GO-KEGG"
#SBATCH --output=log/%x_%A_%a.log
##################################################################

## RNA-seq Workflow - GO-KEGG analysis #############
##
## This workflow will carry out the following steps:
## 1 - GO analysis
## 2 - KEGG analysis

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

## Running workflow

echo "## Starting enrichment analysis ##########"

sed -n ${N}p $FACTORS | while IFS="," read FCT
do
    echo "Processing Factor = ${FCT}"

    Rscript --vanilla code/GO_workflow_whitefly.R ${FCT}
    # Rscript --vanilla code/KEGG_workflow_whitefly_v2.R ${FCT}

done