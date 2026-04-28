#!/bin/bash
#SBATCH --job-name=gopeaks_peak_calling
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=isacaste@iu.edu
#SBATCH --output=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-calling/gopeaks/logs/gopeaks%A_%a.out
#SBATCH --error=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-calling/gopeaks/logs/gopeaks%A_%a.err
#SBATCH --array=0-37%5
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH -A r00750
#SBATCH -p general

# ============================================================
# GoPeaks Peak Calling — CUT&Tag, no control, hg38
# 38 samples: H3K27ac + H3K27me3 (Abbasova et al. 2025 data)
# Yashar et al. 2022 (Genome Biology)
#
# Fixes vs previous version:
#   - Removed --blacklist (not a valid GoPeaks flag)
#   - Mark detection hardcoded from metadata instead of sample name
#   - H3K27ac uses --broad (wide domains, active enhancers)
#   - H3K27me3 uses default narrow settings
#   - Blacklist filtering applied as post-processing step
#
# Input:  data/dedup-bam-files/
# Output: results/peak-calling/gopeaks/<sample>/
# ============================================================

source ~/.bashrc
conda activate gopeaks_env

# --- Paths ---
DEDUP_DIR=/N/project/Krolab/isabella/data/dedup-bam-files
OUT_DIR=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-calling/gopeaks
BLACKLIST=/N/project/Krolab/isabella/data/encode-files/hg38/ENCFF000KJP_hg38.bed

mkdir -p "${OUT_DIR}/logs"

# --- Build sample list ---
mapfile -t SAMPLES < <(ls "${DEDUP_DIR}"/*.bam | xargs -n1 basename | sed 's/\.bam$//')

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
BAM="${DEDUP_DIR}/${SAMPLE}.bam"

echo "============================="
echo "Sample:     ${SAMPLE}"
echo "Input BAM:  ${BAM}"
echo "Array ID:   ${SLURM_ARRAY_TASK_ID}"
echo "Date:       $(date)"
echo "============================="

mkdir -p "${OUT_DIR}/${SAMPLE}"

# --- Histone mark lookup from metadata ---
# H3K27me3 samples (7 total): all others are H3K27ac (31 total)
declare -A MARK_MAP
MARK_MAP["SRR31972717_dedup"]="H3K27me3"
MARK_MAP["SRR31972726_dedup"]="H3K27me3"
MARK_MAP["SRR31972730_dedup"]="H3K27me3"
MARK_MAP["SRR31972735_dedup"]="H3K27me3"
MARK_MAP["SRR31972738_dedup"]="H3K27me3"
MARK_MAP["SRR31972739_dedup"]="H3K27me3"
MARK_MAP["SRR31972748_dedup"]="H3K27me3"

MARK="${MARK_MAP[$SAMPLE]:-H3K27ac}"
echo "Histone mark: ${MARK}"

# --- Set GoPeaks flags based on mark ---
# H3K27ac:  --broad (step=5000, slide=1000 — wide enhancer domains)
# H3K27me3: default narrow settings (step=100, slide=50)
if [[ "${MARK}" == "H3K27ac" ]]; then
    echo "Using --broad mode for H3K27ac"
    GOPEAKS_FLAGS="--broad"
else
    echo "Using default narrow mode for H3K27me3"
    GOPEAKS_FLAGS=""
fi

# --- Ensure BAM is indexed ---
if [ ! -f "${BAM}.bai" ]; then
    echo "Index not found — indexing BAM..."
    samtools index "${BAM}"
fi

# --- Step 1: GoPeaks peak calling ---
gopeaks \
    --bam "${BAM}" \
    --prefix "${OUT_DIR}/${SAMPLE}/${SAMPLE}" \
    ${GOPEAKS_FLAGS}

echo "GoPeaks peak calling done — $(date)"

# --- Step 2: Post-hoc blacklist filtering ---
# GoPeaks has no --blacklist flag, so we filter the output BED manually
PEAKS_BED="${OUT_DIR}/${SAMPLE}/${SAMPLE}_peaks.bed"
FILTERED_BED="${OUT_DIR}/${SAMPLE}/${SAMPLE}_peaks_blacklist_filtered.bed"

if [ -f "${PEAKS_BED}" ]; then
    echo "Applying blacklist filter..."
    bedtools intersect \
        -a "${PEAKS_BED}" \
        -b "${BLACKLIST}" \
        -v \
        > "${FILTERED_BED}"
    BEFORE=$(wc -l < "${PEAKS_BED}")
    AFTER=$(wc -l < "${FILTERED_BED}")
    echo "Peaks before blacklist filter: ${BEFORE}"
    echo "Peaks after blacklist filter:  ${AFTER}"
else
    echo "WARNING: Expected peaks BED not found at ${PEAKS_BED}"
    echo "Check GoPeaks output files in ${OUT_DIR}/${SAMPLE}/"
    ls "${OUT_DIR}/${SAMPLE}/"
fi

echo "GoPeaks done for ${SAMPLE} — $(date)"
