#!/bin/bash
#SBATCH --job-name=chromhmm_enrich
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/chromhmm/logs/chromhmm_%A_%a.out
#SBATCH --error=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/chromhmm/logs/chromhmm_%A_%a.err
#SBATCH --array=0-1
#SBATCH --mail-user=isacaste@iu.edu
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH -A r00750

# ChromHMM State Enrichment Analysis
# Overlap peaks (GoPeaks, MACS2) with K562 chromatin states
# Calculate enrichment for each state

module load bedtools

# --- Paths ---
CHROMHMM_FILE="/N/project/Krolab/isabella/data/chromhmm/E123_15_coreMarks_hg38lift_mnemonics.bed"
PEAK_DIR="/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-calling"
OUT_DIR="/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/chromhmm"
METADATA="/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/metadata.xlsx"

mkdir -p "${OUT_DIR}/logs"
mkdir -p "${OUT_DIR}/overlaps"

# Peak caller options (just GoPeaks and MACS2)
TOOLS=("GoPeaks" "MACS2")
TOOL="${TOOLS[$SLURM_ARRAY_TASK_ID]}"

echo "============================="
echo "Tool:       ${TOOL}"
echo "Array ID:   ${SLURM_ARRAY_TASK_ID}"
echo "============================="

# Check ChromHMM file exists
if [ ! -f "${CHROMHMM_FILE}" ]; then
    echo "ERROR: ChromHMM file not found: ${CHROMHMM_FILE}"
    echo "Run download_chromhmm.sh first"
    exit 1
fi

echo "ChromHMM states: $(wc -l < ${CHROMHMM_FILE}) regions"

# Output file
OUTPUT_CSV="${OUT_DIR}/${TOOL}_chromhmm_enrichment.csv"

# Write CSV header
echo "Sample,Histone_Mark,Tool,State,Peak_Count,State_Bp,Overlap_Bp,Enrichment_Score" > "${OUTPUT_CSV}"

# Representative samples for each histone mark (to keep runtime reasonable)
# H3K27ac: 5 samples (different antibodies)
# H3K27me3: All 7 samples (small set)
H3K27AC_SAMPLES=(SRR31972716 SRR31972718 SRR31972719 SRR31972720 SRR31972727)
H3K27ME3_SAMPLES=(SRR31972717 SRR31972726 SRR31972730 SRR31972735 SRR31972738 SRR31972739 SRR31972748)

ALL_SAMPLES=("${H3K27AC_SAMPLES[@]}" "${H3K27ME3_SAMPLES[@]}")

# Histone mark lookup
declare -A MARK_MAP
MARK_MAP["SRR31972716"]="H3K27ac"
MARK_MAP["SRR31972717"]="H3K27me3"
MARK_MAP["SRR31972718"]="H3K27ac"
MARK_MAP["SRR31972719"]="H3K27ac"
MARK_MAP["SRR31972720"]="H3K27ac"
MARK_MAP["SRR31972726"]="H3K27me3"
MARK_MAP["SRR31972727"]="H3K27ac"
MARK_MAP["SRR31972730"]="H3K27me3"
MARK_MAP["SRR31972735"]="H3K27me3"
MARK_MAP["SRR31972738"]="H3K27me3"
MARK_MAP["SRR31972739"]="H3K27me3"
MARK_MAP["SRR31972748"]="H3K27me3"

# Get total genome size covered by each chromatin state
echo "Calculating chromatin state sizes..."
STATE_SIZES=$(mktemp)
awk '{len=$3-$2; state=$4; bp[state]+=len} END {for (s in bp) print s"\t"bp[s]}' "${CHROMHMM_FILE}" > "${STATE_SIZES}"

# Total genome size in chromatin states
TOTAL_BP=$(awk '{sum+=$2} END {print sum}' "${STATE_SIZES}")
echo "Total annotated genome: ${TOTAL_BP} bp"

# Process each sample
for SAMPLE in "${ALL_SAMPLES[@]}"; do
    
    MARK="${MARK_MAP[$SAMPLE]}"
    echo ""
    echo "Processing ${SAMPLE} (${MARK})..."
    
    # Get peak file based on tool
    if [ "${TOOL}" == "MACS2" ]; then
        PEAK_FILE="${PEAK_DIR}/macs2/${SAMPLE}_dedup/${SAMPLE}_dedup_peaks.narrowPeak"
    elif [ "${TOOL}" == "GoPeaks" ]; then
        PEAK_FILE="${PEAK_DIR}/gopeaks/${SAMPLE}_dedup/${SAMPLE}_dedup_peaks.bed"
    fi
    
    if [ ! -f "${PEAK_FILE}" ]; then
        echo "WARNING: Peak file not found: ${PEAK_FILE}"
        continue
    fi
    
    PEAK_COUNT=$(wc -l < "${PEAK_FILE}")
    echo "  Peaks: ${PEAK_COUNT}"
    
    # Intersect peaks with chromatin states
    OVERLAP_FILE="${OUT_DIR}/overlaps/${SAMPLE}_${TOOL}_chromhmm_overlap.bed"
    
    bedtools intersect \
        -a "${PEAK_FILE}" \
        -b "${CHROMHMM_FILE}" \
        -wa -wb \
        > "${OVERLAP_FILE}"
    
    # Calculate enrichment for each state
    # Enrichment = (observed overlap bp / total peak bp) / (state size / total genome bp)
    
    TOTAL_PEAK_BP=$(awk '{sum+=$3-$2} END {print sum}' "${PEAK_FILE}")
    
    # For each chromatin state, calculate overlap
    while read STATE_NAME STATE_BP; do
        
        # Get overlap bp for this state
        OVERLAP_BP=$(awk -v state="${STATE_NAME}" '$NF==state {sum+=$3-$2} END {print sum+0}' "${OVERLAP_FILE}")
        
        # Calculate enrichment score
        # Enrichment = (overlap_bp / total_peak_bp) / (state_bp / total_bp)
        if [ "${TOTAL_PEAK_BP}" -gt 0 ] && [ "${TOTAL_BP}" -gt 0 ]; then
            ENRICHMENT=$(awk -v overlap="${OVERLAP_BP}" \
                             -v total_peak="${TOTAL_PEAK_BP}" \
                             -v state_bp="${STATE_BP}" \
                             -v total_bp="${TOTAL_BP}" \
                             'BEGIN {
                                 observed = overlap / total_peak;
                                 expected = state_bp / total_bp;
                                 if (expected > 0) {
                                     enrichment = observed / expected;
                                 } else {
                                     enrichment = 0;
                                 }
                                 printf "%.4f", enrichment
                             }')
        else
            ENRICHMENT="0.0000"
        fi
        
        echo "${SAMPLE},${MARK},${TOOL},${STATE_NAME},${PEAK_COUNT},${STATE_BP},${OVERLAP_BP},${ENRICHMENT}" >> "${OUTPUT_CSV}"
        
    done < "${STATE_SIZES}"
    
    echo "  ✓ Enrichment calculated"
    
done

# Cleanup
rm -f "${STATE_SIZES}"

echo "============================="
echo "ChromHMM enrichment complete for ${TOOL}"
echo "Output: ${OUTPUT_CSV}"
echo "============================="
