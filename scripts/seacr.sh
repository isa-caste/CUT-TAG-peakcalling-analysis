#!/bin/bash
#SBATCH --job-name=seacr_peak_calling
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=isacaste@iu.edu
#SBATCH --output=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-calling/seacr/logs/seacr_%A_%a.out
#SBATCH --error=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-calling/seacr/logs/seacr_%A_%a.err
#SBATCH --array=0-37%5
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH -A r00750
#SBATCH -p general

# ============================================================
# SEACR Peak Calling — CUT&Tag, no control, hg38
# 38 samples: H3K27ac + H3K27me3 (Abbasova et al. 2025 data)
#
# Matches paper exactly:
#   SEACR_1.3.sh input_bedgraph 0.01 non stringent out_name
#
# Input:  data/bed-files/
# Output: results/peak-calling/seacr/<sample>/
# ============================================================

source ~/.bashrc
conda activate seacr_env

# --- Paths ---
SEACR=/N/project/Krolab/isabella/tools/SEACR/SEACR_1.3.sh
BED_DIR=/N/project/Krolab/isabella/data/bed-files
OUT_DIR=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-calling/seacr
GENOME_SIZES=/N/project/Krolab/isabella/H3K9me2-Research/annotations/hg38_primary_chrom_sizes.txt

mkdir -p "${OUT_DIR}/logs"

# --- Build sample list ---
mapfile -t SAMPLES < <(ls "${BED_DIR}"/*.bed | xargs -n1 basename | sed 's/\.bed$//')

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
BED="${BED_DIR}/${SAMPLE}.bed"

echo "============================="
echo "Sample:     ${SAMPLE}"
echo "Input BED:  ${BED}"
echo "Array ID:   ${SLURM_ARRAY_TASK_ID}"
echo "Date:       $(date)"
echo "============================="

mkdir -p "${OUT_DIR}/${SAMPLE}"

# --- Step 1: Generate fragment bedgraph (skip if already exists) ---
BEDGRAPH="${OUT_DIR}/${SAMPLE}/${SAMPLE}.fragments.bedgraph"

if [ -s "${BEDGRAPH}" ]; then
    echo "Bedgraph already exists and is non-empty — skipping generation"
else
    echo "Generating bedgraph for ${SAMPLE}..."
    # Filter to primary chromosomes only (chr1-22, chrX, chrY, chrM)
    # Unplaced contigs (GL*, KI*, etc.) are not in hg38_primary_chrom_sizes.txt
    # and will cause bedtools genomecov to exit with an error
    grep -E '^chr([0-9]+|X|Y|M)\s' "${BED}" | \
        sort -k1,1 -k2,2n | \
        bedtools genomecov \
            -bg \
            -i stdin \
            -g "${GENOME_SIZES}" \
        > "${BEDGRAPH}"
    echo "Bedgraph done — $(date)"
fi

# --- Step 2: SEACR peak calling (matching Abbasova et al. 2025) ---
# Argument order: <bedgraph> <threshold> <norm> <mode> <output_prefix>
#   0.01       — top 1% signal threshold (no-control mode)
#   non        — no normalisation (no IgG control available)
#   stringent  — peak calling mode (intersect of two criteria)
bash "${SEACR}" \
    "${BEDGRAPH}" \
    0.01 \
    non \
    stringent \
    "${OUT_DIR}/${SAMPLE}/${SAMPLE}"

echo "SEACR done for ${SAMPLE} — $(date)"

conda deactivate
