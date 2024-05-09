#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00     # 4 hours
#SBATCH --job-name="conda"
#SBATCH --output=log/%x_%j.log
#SBATCH -p short # This is the default partition, you can use any of the following; epyc, short, intel, batch, highmem, gpu

conda create --name DEG-analysis --file code/conda_spec_file.txt