#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00     # 6 hours
#SBATCH --job-name="trim_galore"
#SBATCH --output=log/%x_%j.log
#SBATCH -p epyc # This is the default partition, you can use any of the following; intel, batch, highmem, gpu

################### Trim_galore.sh ###################

# Load modules
module load trim_galore
module load fastqc


# Reading in variables
SAMPLENAME=$1
FQ1=$2
FQ2=$3
OUT_DIR=analysis/trim_galore

# Create analysis folder (if it doesn't exist)
[ ! -d "$OUT_DIR" ] && mkdir -p "$OUT_DIR"

# Running trim_galore
trim_galore --fastqc --gzip --trim-n -j 4 --paired -o $OUT_DIR --basename ${SAMPLENAME} ${FQ1} ${FQ2}

# Running trim_galore (single-end)
# trim_galore --fastqc --gzip --trim-n -j 4 -o $OUT_DIR --basename ${SAMPLENAME} ${FQ1}
