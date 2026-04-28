#!/bin/bash
#SBATCH --job-name=seacr_nonstringent
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/seacr_nonstringent_%A_%a.out
#SBATCH --error=logs/seacr_nonstringent_%A_%a.err
#SBATCH --array=0-37
#SBATCH --mail-user=isacaste@iu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A r00750

source ~/.bashrc
module load conda
conda activate seacr_env

# --- Paths ---
SEACR=/N/project/Krolab/isabella/tools/SEACR/SEACR_1.3.sh
BED_DIR=/N/project/Krolab/isabella/data/bed-files
OUT_DIR=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-calling/seacr
GENOME_SIZES=/N/project/Krolab/isabella/H3K9me2-Research/annotations/hg38_primary_chrom_sizes.txt

mkdir -p "${OUT_DIR}/logs"
 
# --- Build sample list ---
mapfile -t SAMPLES < <(ls -d "${OUT_DIR}"/*_fragments | xargs -n1 basename)
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
 
echo "============================="
echo "Sample:     ${SAMPLE}"
echo "Input BED:  ${BED}"
echo "Array ID:   ${SLURM_ARRAY_TASK_ID}"
echo "============================="
 
mkdir -p "${OUT_DIR}/${SAMPLE}"
 
# --- Step 1: Check if bedgraph exists ---
BEDGRAPH="${OUT_DIR}/${SAMPLE}/${SAMPLE}.fragments.bedgraph"
 
if [ ! -s "${BEDGRAPH}" ]; then
    echo "ERROR: Bedgraph not found: ${BEDGRAPH}"
    echo "You need to generate bedgraphs first (should exist from stringent run)"
    exit 1
fi
 
echo "Using existing bedgraph: ${BEDGRAPH}"
# --- Step 2: SEACR peak calling - NON-STRINGENT mode ---
# Argument order: <bedgraph> <threshold> <norm> <mode> <output_prefix>
#   0.01          — top 1% signal threshold (no-control mode)
#   non           — no normalisation (no IgG control available)
#   non-stringent — relaxed peak calling mode (union of two criteria)
 
bash "${SEACR}" \
    "${BEDGRAPH}" \
    0.01 \
    non \
    relaxed \
    "${OUT_DIR}/${SAMPLE}/${SAMPLE}"
 
echo "SEACR non-stringent done for ${SAMPLE} — $(date)"

conda deactivate
