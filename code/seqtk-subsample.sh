#!/bin/bash -l 
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=6:00:00
#SBATCH --job-name="seqtk"
#SBATCH --output=log/%x_%j.log
#####################################
# Subsample
######################################


# Load module
module load seqtk

# Set variables
SRR_ID=$1
IN_DIR=raw
SEED=2023

seqtk sample -s ${SEED} $IN_DIR/${SRR_ID}_1.fastq.gz 1000000 > ${IN_DIR}/${SRR_ID}_sub1M_1.fq.gz
seqtk sample -s ${SEED} $IN_DIR/${SRR_ID}_2.fastq.gz 1000000 > ${IN_DIR}/${SRR_ID}_sub1M_2.fq.gz

