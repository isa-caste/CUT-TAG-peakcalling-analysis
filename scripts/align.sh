#!/bin/bash
#SBATCH --mail-user=isacaste@iu.edu
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --job-name=align
#SBATCH -o align.out
#SBATCH -e align.err
#SBATCH --nodes=2
#SBATCH --mem=16g
#SBATCH -p gpu
#SBATCH --ntasks-per-node=2
#SBATCH --gpus-per-node=1
#SBATCH --time=1-23:59:00
#SBATCH -A r00750

# Load required modules
module load bowtie/2.5.1
module load samtools

# define directories
TRIMMED_DIR=/N/project/Krolab/isabella/data/trimmed-files
INDEX_DIR=/N/project/Krolab/isabella/H3K9me2-Research/rna-seq/hg38_bt2_index
GENOME_PREFIX="GRCh38"
OUTDIR=/N/project/Krolab/isabella/data/bam-files

# Move to trimmed fastq directory
cd "$TRIMMED_DIR"

# make output directory it if doesn't exist
mkdir -p "$OUTDIR"

# Run Bowtie2 alignment for paired end reads
for R1_FILE in *_1_val_1.fq.gz; do
    # Check if file exists
    if [ ! -f "$R1_FILE" ]; then
        echo "No R1 files found. Exiting."
        exit 1
    fi

    # Define paired-end filename
    R2_FILE="${R1_FILE/_1_val_1.fq.gz/_2_val_2.fq.gz}"
    BASENAME="${R1_FILE/_1_val_1.fq.gz/}"

    # Verify R2 file exists
    if [ ! -f "$R2_FILE" ]; then
        echo "Warning: $R2_FILE not found. Skipping $BASENAME"
        continue
    fi
    echo "Begining aligment for $BASENAME..."
    # Run Bowtie2 alignment for paired-end reads
    bowtie2 \
        -x "${INDEX_DIR}/${GENOME_PREFIX}" \
        -1 "$R1_FILE" \
        -2 "$R2_FILE" \
        -S "${BASENAME}.sam" \
        --no-mixed \
        --no-discordant \
        -p 4 \
        2> "${BASENAME}_align.err"

    # Convert SAM to sorted BAM and index
    echo "Converting SAM to BAM and sorting..."
    samtools view -bS "${BASENAME}.sam" | samtools sort -m 4G -o "${OUTDIR}/${BASENAME}_sorted.bam"

    # Index the BAM file
    samtools index "${OUTDIR}/${BASENAME}_sorted.bam"

    #  Calculate alignment statistics
    samtools flagstat "${OUTDIR}/${BASENAME}_sorted.bam" > "${OUTDIR}/${BASENAME}_flagstat.txt"

    # Cleanup
    rm "${BASENAME}.sam"

    echo "$BASENAME alignment complete!"
done
