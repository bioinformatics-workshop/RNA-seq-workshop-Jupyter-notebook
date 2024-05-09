#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00     # 4 hours
#SBATCH --job-name="fastqc"
#SBATCH --output=log/%x_%j.log
#SBATCH -p epyc # This is the default partition, you can use any of the following; epyc, short, intel, batch, highmem, gpu

module load fastqc

# Load variables
SRR_ID=$1
IN_DIR=raw
OUT_DIR=analysis/fastqc

# Create analysis folder (if it doesn't exist)
[ ! -d "$OUT_DIR" ] && mkdir -p "$OUT_DIR"

# -t : threads
# --noextract : do not extract zip output files

fastqc -f fastq -t 4 --noextract -o $OUT_DIR $IN_DIR/${SRR_ID}*.fastq.gz

