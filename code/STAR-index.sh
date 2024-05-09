#!/bin/bash -l
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=32
#SBATCH --mem=100g
#SBATCH --time=1-23:00:00
#SBATCH --job-name="STAR-index"
#SBATCH --output=log/%x_%j.log

## Load environment
module load star

FASTA=genome/TAIR10_Chr.all.fasta
GTF=genome/Araport11_GFF3_genes_transposons.201606.gtf
GFF=genome/Araport11_GFF3_genes_transposons.201606.gff
IDX_DIR=index

PARAMS="--runThreadN 32
	--runMode genomeGenerate
        --genomeDir $IDX_DIR
        --genomeFastaFiles $FASTA
        --sjdbGTFfile $GTF
        --sjdbOverhang 99
        --genomeSAindexNbases 12"

STAR $PARAMS
