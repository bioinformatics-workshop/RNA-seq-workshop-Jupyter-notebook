#!/bin/bash -l
#SBATCH --partition=epyc
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --mem=100g
#SBATCH --time=1-12:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name="01-rnaseq-align"
#SBATCH --output=log/%x_%A_%a.log
##################################################################

## RNA-seq Workflow - Alignment #############
##
## This workflow will carry out the following steps:
## 1 - Run QC
## 2 - Trim adapters and low-quality bases
## 3 - Align reads using STAR
## 4 - Generate counts

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
LIST_DIR=( ${FASTQC_DIR}
           ${TRIMGALORE_DIR}
           ${STAR_DIR}/${SAMPLENAME} 
           ${FEATURECOUNTS_DIR} )

for DIR in "${LIST_DIR[@]}"
do
    [ ! -d "$DIR" ] && mkdir -p "${DIR}"
done

## Running workflow

echo "## Starting RNA-seq Workflow ##########"

sed -n ${N}p $SAMPLES | while IFS="," read SAMPLENAME FASTQ1 FASTQ2
do
    echo "Processing sample = ${SAMPLENAME}"
    echo "Paired-end = ${PAIRED}"
    echo "Strandedness = ${STRANDED}"

    run_fastqc
    run_trimgalore
    run_star
    run_featurecounts
    sam_index
    bam2bw

done

