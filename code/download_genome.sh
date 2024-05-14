#!/bin/bash
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=6:00:00
#SBATCH --job-name="genome-download"
#SBATCH --output=log/%x_%j.log
#####################################
# Download reference genomes
######################################


# Arabidopsis reference genome (TAIR10/Araport11)
# Resources available @ https://www.arabidopsis.org

# Set Parameters
CHROM_URL=https://www.arabidopsis.org/download/file?path=Sequences/Assemblies/TAIR9_chr_all.fas
GFF_URL=https://www.arabidopsis.org/download/file?path=Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.current.gff.gz
GTF_URL=https://www.arabidopsis.org/download/file?path=Genes/Araport11_genome_release/Araport11_GTF_genes_transposons.current.gtf.gz
README_URL=https://www.arabidopsis.org/download/file?path=Genes/Araport11_genome_release/README.202404.md
OUT_DIR=genome  


# Download each link
for URL in CHROM_URL GFF_URL GTF_URL README_URL
do
    wget -r --reject="index.html*" --no-parent -l1 -nH --directory-prefix=${OUT_DIR} ${!URL}
done




