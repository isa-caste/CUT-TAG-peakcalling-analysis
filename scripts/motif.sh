#!/bin/bash
#SBATCH --job-name=motif_analysis
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/$
#SBATCH --error=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/m$
#SBATCH --array=0-3
#SBATCH --mail-user=isacaste@iu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A r00750

# HOMER Motif Analysis- find enriched transcription factor binding motifs in peaks

module load bedtools

# Paths
PEAK_DIR="/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-calling"
OUT_DIR="/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/motif-analysis"
GENOME="/N/project/Krolab/isabella/data/genomes/hg38/hg38.fa"

mkdir -p "${OUT_DIR}/logs"

# Conditions: Tool_Mark
CONDITIONS=("GoPeaks_H3K27ac" "GoPeaks_H3K27me3" "MACS2_H3K27ac" "MACS2_H3K27me3")
CONDITION="${CONDITIONS[$SLURM_ARRAY_TASK_ID]}"

IFS='_' read -r TOOL MARK <<< "$CONDITION"

echo "============================="
echo "Condition:  ${CONDITION}"
echo "Tool:       ${TOOL}"
echo "Mark:       ${MARK}"
echo "Array ID:   ${SLURM_ARRAY_TASK_ID}"
echo "============================="

# Representative samples (same as ChromHMM/GO)
if [ "${MARK}" == "H3K27ac" ]; then
    SAMPLES=(SRR31972716 SRR31972718 SRR31972719 SRR31972720 SRR31972727)
else
    SAMPLES=(SRR31972717 SRR31972726 SRR31972730 SRR31972735 SRR31972738 SRR31972739 SRR31972748)
fi

# Merge peaks from all representative samples
MERGED_PEAKS="${OUT_DIR}/${CONDITION}_merged_peaks.bed"
rm -f "${MERGED_PEAKS}"

for SAMPLE in "${SAMPLES[@]}"; do
    if [ "${TOOL}" == "MACS2" ]; then
        PEAK_FILE="${PEAK_DIR}/macs2/${SAMPLE}_dedup/${SAMPLE}_dedup_peaks.narrowPeak"
    else
        PEAK_FILE="${PEAK_DIR}/gopeaks/${SAMPLE}_dedup/${SAMPLE}_dedup_peaks.bed"
    fi
    
    if [ -f "${PEAK_FILE}" ]; then
        # Extract chr, start, end only and filter to standard chromosomes
        cut -f1-3 "${PEAK_FILE}" | grep -E '^chr[0-9XY]+\s' >> "${MERGED_PEAKS}"
    fi
done

# Sort and merge overlapping peaks
MERGED_SORTED="${OUT_DIR}/${CONDITION}_merged_sorted.bed"
sort -k1,1 -k2,2n "${MERGED_PEAKS}" | bedtools merge -i stdin > "${MERGED_SORTED}"

PEAK_COUNT=$(wc -l < "${MERGED_SORTED}")
echo "Merged peak regions: ${PEAK_COUNT}"

if [ "${PEAK_COUNT}" -lt 100 ]; then
    echo "WARNING: Too few peaks for motif analysis (< 100). Skipping."
    exit 0
fi

# Check if HOMER is available
if ! command -v findMotifsGenome.pl &> /dev/null; then
    echo "ERROR: HOMER not found. Trying to load module..."
    module load homer 2>/dev/null || true
    
    if ! command -v findMotifsGenome.pl &> /dev/null; then
        echo "ERROR: HOMER still not available. Install or load module."
        exit 1
    fi
fi

# Check if genome exists
if [ ! -f "${GENOME}" ]; then
    echo "WARNING: Genome file not found: ${GENOME}"
    echo "HOMER will download hg38 automatically"
    GENOME="hg38"
fi

# Run HOMER motif analysis
MOTIF_OUT="${OUT_DIR}/${CONDITION}_motifs"
mkdir -p "${MOTIF_OUT}"
PREPARSE_DIR="${OUT_DIR}/preparsed"
mkdir -p "${PREPARSE_DIR}"

echo ""
echo "Running HOMER motif analysis..."

findMotifsGenome.pl \
    "${MERGED_SORTED}" \
    ${GENOME} \
    "${MOTIF_OUT}" \
    -size 200 \
    -mask \
    -p 8 \
    -preparsedDir "${PREPARSE_DIR}"

echo "Motif analysis complete for ${CONDITION}"
echo "Output: ${MOTIF_OUT}"

# Extract top motifs summary
if [ -f "${MOTIF_OUT}/knownResults.txt" ]; then
    echo ""
    echo "Top 10 Known Motifs:"
    head -11 "${MOTIF_OUT}/knownResults.txt" | tail -10 | \
        awk 'BEGIN{OFS="\t"} {print $1, $2, $6}'
fi
