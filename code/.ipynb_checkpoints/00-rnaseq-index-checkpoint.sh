#!/bin/bash -l
#SBATCH --partition=epyc
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --mem=100g
#SBATCH --time=1-12:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name="00-rnaseq-index"
#SBATCH --output=log/%x_%j.log
##################################################################

## RNA-seq Workflow - Genome Indexing #############

## Get configuration file
if [ -f code/config.txt ]; then
  source code/config.txt
else
    echo "You're missing the config.txt file"
    exit
fi

## Load modules and environments
module load ${MODULE_LIST}

NUM_CORES=8

STAR --runThreadN ${NUM_CORES} \
     --runMode genomeGenerate \
     --genomeDir ${INDEX} \
     --genomeFastaFiles ${GENOME_FASTA} \
     --sjdbGTFfile $GTF \
     --sjdbOverhang 99 \
     --genomeSAindexNbases 12