#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00     # 4 hours
#SBATCH --job-name="sra-download"
#SBATCH -p epyc # This is the default partition, you can use any of the following; epyc, short, intel, batch, highmem, gpu
#SBATCH --output=log/%x_%j.log

#################################

# current version 3.0.0 4/28/22
module load sratoolkit

# Get parameters
SRR_ID=$1
OUT_DIR=raw

fastq-dump --split-3 --gzip --outdir ${OUT_DIR} ${SRR_ID}
