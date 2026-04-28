#!/bin/bash
#SBATCH --mail-user=isacaste@iu.edu
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --nodes=2
#SBATCH --mem=16g
#SBATCH -p gpu
#SBATCH --ntasks-per-node=2
#SBATCH --gpus-per-node=1
#SBATCH --time=1-23:59:00
#SBATCH --job-name=trimming
#SBATCH -o trim.out
#SBATCH -e trim.err
#SBATCH -A r00750

# load conda env
module load conda
conda activate trimming

# set directories
INPUT_DIR=/N/project/Krolab/isabella/data/fastq-files
OUTPUT_DIR=/N/project/Krolab/isabella/data/trimmed-files
MQC_DIR=/N/project/Krolab/isabella/data/trim-mqc

mkdir -p $OUTPUT_DIR $MQC_DIR

cd $INPUT_DIR

for r1 in *_1.fastq.gz; do
    r2="${r1/_1.fastq.gz/_2.fastq.gz}"
    trim_galore --paired --fastqc --nextera --length 20 -o $OUTPUT_DIR $r1 $r2
done

# move fastqc files to mqc directory
mv $OUTPUT_DIR/*_fastqc.html $OUTPUT_DIR/*_fastqc.zip $MQC_DIR/

# run multiqc
multiqc $MQC_DIR -o $MQC_DIR

# deactivate conda env
conda deactivate

