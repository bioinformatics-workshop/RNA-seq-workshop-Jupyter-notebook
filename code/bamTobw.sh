#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --time=02:00:00
#SBATCH --job-name="bam2bw"
#SBATCH --output=log/%x_%j.log
#SBATCH -p epyc # This is the default partition, you can use any of the following; intel, batch, highmem, gpu

# Shell script to generate BigWig files from STAR bam output
# Since the dataset is non-stranded, we can convert directly from the BAM file

# Load deeptools
module load deeptools

# Set variables
SAMPLENAME=$1
IN_DIR=analysis/star

bamCoverage -b $IN_DIR/${SAMPLENAME}_Aligned.sortedByCoord.out.bam \
    -o $IN_DIR/${SAMPLENAME}.bw \
    -p 4 \
    --normalizeUsing CPM \
    --exactScaling \
    --skipNAs
