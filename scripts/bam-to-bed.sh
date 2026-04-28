#!/bin/bash
#SBATCH --job-name=bam_to_bed
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=isacaste@iu.edu
#SBATCH --output=bamtobed_%A_%a.out
#SBATCH --error=bamtobed_%A_%a.err
#SBATCH --array=0-37%5
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH -A r00750
#SBATCH -p general

source ~/.bashrc
module load bamtools
module load samtools

# --- Paths ---
BAM_DIR=/N/project/Krolab/isabella/data/dedup-bam-files
BED_DIR=/N/project/Krolab/isabella/data/bed-files
TMP_DIR=/tmp/namesorted_bams

mkdir -p "${BED_DIR}" "${TMP_DIR}"

# --- Build sample list ---
mapfile -t SAMPLES < <(ls "${BAM_DIR}"/*_dedup.bam | xargs -n1 basename | sed 's/_dedup\.bam$//')

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
BAM="${BAM_DIR}/${SAMPLE}_dedup.bam"
BED="${BED_DIR}/${SAMPLE}_fragments.bed"
NAMESORTED="${TMP_DIR}/${SAMPLE}_namesorted.bam"

echo "Sample:     ${SAMPLE}"
echo "Input BAM:  ${BAM}"
echo "Output BED: ${BED}"

# --- Step 1: Name-sort BAM (required for -bedpe mate pairing) ---
echo "Name-sorting BAM..."
samtools sort -n \
    -@ "${SLURM_CPUS_PER_TASK}" \
    "${BAM}" \
    -o "${NAMESORTED}"

echo "Name-sort done — $(date)"

# --- Step 2: Convert to fragment BED ---
# -bedpe: outputs one line per fragment (both mates)
# awk:    extracts chr, fragment_start, fragment_end (cols 1, 2, 6)
#         and filters out unmapped reads (. -1 -1)
echo "Converting to fragment BED..."
bedtools bamtobed \
    -i "${NAMESORTED}" \
    -bedpe 2>/dev/null | \
    awk 'BEGIN {OFS="\t"} $1 != "." && $2 >= 0 && $6 >= 0 {print $1, $2, $6}' | \
    sort -k1,1 -k2,2n \
    > "${BED}"

FRAG_COUNT=$(wc -l < "${BED}")
echo "Fragment BED done — ${FRAG_COUNT} fragments — $(date)"

# --- Step 3: Clean up name-sorted BAM to save space ---
rm -f "${NAMESORTED}"
echo "Cleaned up tmp file"

echo "Done: ${SAMPLE}"
