#!/bin/bash
#SBATCH --mail-user=isacaste@iu.edu
#SBATCH --nodes=2
#SBATCH --mem=16g
#SBATCH -p gpu
#SBATCH --ntasks-per-node=2
#SBATCH --gpus-per-node=1
#SBATCH --time=1-23:59:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --job-name=qc
#SBATCH -o qc.out
#SBATCH -e qc.err
#SBATCH -A r00750

# load conda environment
module load conda 
conda activate align-qc-env

# Define input and output directories
INPUT_DIR=/N/project/Krolab/isabella/data/fastq-files
OUTPUT_DIR=/N/project/Krolab/isabella/data/multiqc-fastqc

# Create output dir if it doesn't exist
mkdir -p $OUTPUT_DIR

cd $INPUT_DIR

for file in *.fastq.gz; do
    fastqc -o $OUTPUT_DIR "$file"
done

cd $OUTPUT_DIR
multiqc .
