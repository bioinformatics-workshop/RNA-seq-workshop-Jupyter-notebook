#!/bin/bash -l
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=16
#SBATCH --mem=100g
#SBATCH --time=0-8:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name="STAR-align"
#SBATCH --output=log/%x_%j.log
##################################################################

## Load STAR environment
module load star
module load samtools

## Setting variables
SAMPLENAME=$1
IN_DIR=analysis/trim_galore
OUT_DIR=analysis/star
IDX_DIR=index

# Create analysis folder (if it doesn't exist)
[ ! -d "$OUT_DIR" ] && mkdir -p "$OUT_DIR"

## Setting STAR parameters

PARAMS=" --runThreadN 16
         --runMode alignReads
         --genomeDir ${IDX_DIR}
         --readFilesCommand zcat
         --outSAMtype BAM SortedByCoordinate
         --outFileNamePrefix ${OUT_DIR}/${SAMPLENAME}_
         --outFilterMismatchNmax 0
         --outFilterMultimapNmax 1
         --quantMode GeneCounts
         --outWigType wiggle"

## Running STAR
STAR $PARAMS --readFilesIn $IN_DIR/${SAMPLENAME}_val_1.fq.gz $IN_DIR/${SAMPLENAME}_val_2.fq.gz

## Indexing bam files
samtools index $OUT_DIR/${SAMPLENAME}_*bam
