#!/bin/bash
#SBATCH --job-name=macs2_peak_calling
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=isacaste@iu.edu
#SBATCH --output=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-calling/macs2/logs/macs2%A_%a.out
#SBATCH --error=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-calling/macs2/logs/macs2%A_%a.err
#SBATCH --array=1-37%5
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH -A r00750
#SBATCH -p general

module load conda
source ~/.bashrc
conda activate macs2_env

# --- Paths ---
DEDUP_DIR=/N/project/Krolab/isabella/data/dedup-bam-files
OUT_DIR=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-calling/macs2

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

# --- MACS2 peak calling (matching Abbasova et al. 2025) ---
macs2 callpeak \
    -t "${BAM}" \
    -n "${SAMPLE}" \
    -f BAMPE \
    -g hs \
    -q 1e-5 \
    --keep-dup all \
    --nolambda \
    --nomodel \
    --outdir "${OUT_DIR}/${SAMPLE}"

echo "MACS2 done for ${SAMPLE} — $(date)"
