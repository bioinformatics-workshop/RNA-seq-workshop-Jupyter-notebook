#!/bin/bash
#SBATCH --partition=batch
#SBATCH --cpus-per-task=2
#SBATCH --mem=12g
#SBATCH --time=0-4:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name="featurecount"
#SBATCH --output=log/%x_%j.log

##################################################################
## Shell script to quantify read counts with FeatureCounts

## Load modules
module load subread # v2.0.3

## Setting variables
IN_DIR=analysis/star
OUT_DIR=analysis/featurecounts
STRAND=0 # (0-unstrand, 1-strand, 2-reverse strand)
GTF=genome/Araport11_GFF3_genes_transposons.201606.gtf
FULL_MATRIX_FILE=/rhome/brandonle/bigdata/workshops/demo_files/analysis/featurecounts/featurecounts.txt

# Create analysis folder (if it doesn't exist)
[ ! -d "$OUT_DIR" ] && mkdir -p "$OUT_DIR"

# Create a link to the full featurecounts table
ln -s ${FULL_MATRIX_FILE} ${OUT_DIR}/featurecounts_full.txt

## Set featuercounts params
# -p : indicates paired-end reads
# --countReadPairs : count read pairs instead of reads
# --primary : count primary alignments only

PARAMS="-T 2 \
        -s $STRAND \
        -a $GTF \
        -p \
        -t exon \
        -g gene_id \
        --countReadPairs \
        --primary \
        -o ${OUT_DIR}/featurecounts.txt"

## Running featureCounts
featureCounts ${PARAMS} $IN_DIR/*bam
