#!/bin/bash
#SBATCH --job-name=gopeaks_h3k27me3_fix
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-calling/gopeaks/logs/gopeaks_fix_%A_%a.out
#SBATCH --error=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-calling/gopeaks/logs/gopeaks_fix_%A_%a.err
#SBATCH --array=0-6
#SBATCH --mail-user=isacaste@iu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A r00750

# ============================================================
# GoPeaks Peak Calling — FIX for H3K27me3 samples
# Re-run the 7 failed H3K27me3 samples with --broad flag
# (GoPeaks narrow mode has a bug, but broad mode works)
# ============================================================

source ~/.bashrc
conda activate gopeaks_env

# --- Paths ---
DEDUP_DIR=/N/project/Krolab/isabella/data/dedup-bam-files
OUT_DIR=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-calling/gopeaks
BLACKLIST=/N/project/Krolab/isabella/data/encode-files/hg38/ENCFF000KJP_hg38.bed

mkdir -p "${OUT_DIR}/logs"

# --- Only the 7 H3K27me3 samples that failed ---
H3K27ME3_SAMPLES=(
    SRR31972717_dedup
    SRR31972726_dedup
    SRR31972730_dedup
    SRR31972735_dedup
    SRR31972738_dedup
    SRR31972739_dedup
    SRR31972748_dedup
)

SAMPLE="${H3K27ME3_SAMPLES[$SLURM_ARRAY_TASK_ID]}"
BAM="${DEDUP_DIR}/${SAMPLE}.bam"

echo "============================="
echo "Sample:     ${SAMPLE}"
echo "Input BAM:  ${BAM}"
echo "Array ID:   ${SLURM_ARRAY_TASK_ID}"
echo "Date:       $(date)"
echo "============================="

mkdir -p "${OUT_DIR}/${SAMPLE}"

# --- Ensure BAM is indexed ---
if [ ! -f "${BAM}.bai" ]; then
    echo "Index not found — indexing BAM..."
    samtools index "${BAM}"
fi

# --- Step 1: GoPeaks peak calling with --broad flag ---
# Using --broad for H3K27me3 (narrow mode crashes due to GoPeaks bug)
echo "Running GoPeaks with --broad flag for H3K27me3..."

gopeaks \
    --bam "${BAM}" \
    --prefix "${OUT_DIR}/${SAMPLE}/${SAMPLE}" \
    --broad

if [ $? -ne 0 ]; then
    echo "ERROR: GoPeaks failed for ${SAMPLE}"
    exit 1
fi

echo "GoPeaks peak calling done — $(date)"

# --- Step 2: Post-hoc blacklist filtering ---
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
    echo "ERROR: Expected peaks BED not found at ${PEAKS_BED}"
    ls "${OUT_DIR}/${SAMPLE}/"
    exit 1
fi


