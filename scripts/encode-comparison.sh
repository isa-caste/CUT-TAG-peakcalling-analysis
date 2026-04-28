#!/bin/bash
#SBATCH --job-name=encode_compare
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-comparison/logs/encode_%A_%a.out
#SBATCH --error=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-comparison/logs/encode_%A_%a.err
#SBATCH --array=0-37
#SBATCH --mail-user=isacaste@iu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A r00750

# ENCODE Peak Recovery Analysis
# Compares each peak caller's output against ENCODE gold standard peaks
# Calculates: overlap count, recall %, precision %

module load bedtools

# set up paths
MACS2_DIR=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-calling/macs2
GOPEAKS_DIR=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-calling/gopeaks
SEACR_DIR=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-calling/seacr
ENCODE_DIR=/N/project/Krolab/isabella/data/encode-files/hg38
OUT_DIR=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-comparison/encode-overlap
METADATA=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/metadata.xlsx
 
mkdir -p "${OUT_DIR}/logs"
mkdir -p "${OUT_DIR}/overlaps"
 
# ENCODE reference files
ENCODE_H3K27AC="${ENCODE_DIR}/ENCFF044JNJ_hg38.bed"
ENCODE_H3K27ME3="${ENCODE_DIR}/ENCFF000BXB_hg38.bed"

# Sample list
samples=(SRR31972716 SRR31972717 SRR31972718 SRR31972719 SRR31972720 \
         SRR31972721 SRR31972722 SRR31972723 SRR31972724 SRR31972725 \
         SRR31972726 SRR31972727 SRR31972728 SRR31972729 SRR31972730 \
         SRR31972731 SRR31972732 SRR31972733 SRR31972734 SRR31972735 \
         SRR31972736 SRR31972737 SRR31972738 SRR31972739 SRR31972740 \
         SRR31972741 SRR31972742 SRR31972743 SRR31972744 SRR31972745 \
         SRR31972746 SRR31972747 SRR31972748 SRR31972749 SRR31972750 \
         SRR31972751 SRR31972752 SRR31972753)
# For now, use a simple lookup based on sample ID
declare -A histone_map
histone_map=(
    [SRR31972716]="H3K27ac"
    [SRR31972717]="H3K27me3"
    [SRR31972718]="H3K27ac"
    [SRR31972719]="H3K27ac"
    [SRR31972720]="H3K27ac"
    [SRR31972721]="H3K27ac"
    [SRR31972722]="H3K27ac"
    [SRR31972723]="H3K27ac"
    [SRR31972724]="H3K27ac"
    [SRR31972725]="H3K27ac"
    [SRR31972726]="H3K27me3"
    [SRR31972727]="H3K27ac"
    [SRR31972728]="H3K27ac"
    [SRR31972729]="H3K27ac"
    [SRR31972730]="H3K27me3"
    [SRR31972731]="H3K27ac"
    [SRR31972732]="H3K27ac"
    [SRR31972733]="H3K27ac"
    [SRR31972734]="H3K27ac"
    [SRR31972735]="H3K27me3"
    [SRR31972736]="H3K27ac"
    [SRR31972737]="H3K27ac"
    [SRR31972738]="H3K27me3"
    [SRR31972739]="H3K27me3"
    [SRR31972740]="H3K27ac"
    [SRR31972741]="H3K27ac"
    [SRR31972742]="H3K27ac"
    [SRR31972743]="H3K27ac"
    [SRR31972744]="H3K27ac"
    [SRR31972745]="H3K27ac"
    [SRR31972746]="H3K27ac"
    [SRR31972747]="H3K27ac"
    [SRR31972748]="H3K27me3"
    [SRR31972749]="H3K27ac"
    [SRR31972750]="H3K27ac"
    [SRR31972751]="H3K27ac"
    [SRR31972752]="H3K27ac"
    [SRR31972753]="H3K27ac"
)
SAMPLE="${samples[$SLURM_ARRAY_TASK_ID]}"
HISTONE="${histone_map[$SAMPLE]}"
 
# Select appropriate ENCODE reference
if [ "${HISTONE}" == "H3K27ac" ]; then
    ENCODE_REF="${ENCODE_H3K27AC}"
else
    ENCODE_REF="${ENCODE_H3K27ME3}"
fi
# Get ENCODE peak count
ENCODE_COUNT=$(wc -l < "${ENCODE_REF}")
echo "ENCODE peaks: ${ENCODE_COUNT}"
 
# Function to calculate overlap metrics
calculate_overlap() {
    local tool=$1
    local peaks=$2
    local output_prefix="${OUT_DIR}/overlaps/${SAMPLE}_${tool}"
    
    if [ ! -f "${peaks}" ]; then
        echo "WARNING: Peak file not found: ${peaks}"
        echo "${SAMPLE},${HISTONE},${tool},0,0,0.00,0.00" >> "${OUT_DIR}/encode_comparison_${HISTONE}.csv"
        return
    fi
    
    # Count peaks in query
    query_count=$(wc -l < "${peaks}")
    
    # Find overlaps with ENCODE (at least 1bp overlap)
    bedtools intersect -a "${peaks}" -b "${ENCODE_REF}" -u > "${output_prefix}_overlapping_peaks.bed"
    overlap_count=$(wc -l < "${output_prefix}_overlapping_peaks.bed")
    # Calculate recall: what % of ENCODE peaks did we recover?
    recall=$(awk "BEGIN {printf \"%.2f\", (${overlap_count}/${ENCODE_COUNT})*100}")
    
    # Calculate precision: what % of called peaks overlap ENCODE?
    if [ "${query_count}" -gt 0 ]; then
        precision=$(awk "BEGIN {printf \"%.2f\", (${overlap_count}/${query_count})*100}")
    else
        precision="0.00"
    fi
    
echo "${tool}: ${overlap_count}/${query_count} peaks overlap ENCODE (Recall: ${recall}%, Precision: ${precision}%)"
    
    # Save to CSV
    echo "${SAMPLE},${HISTONE},${tool},${query_count},${overlap_count},${recall},${precision}" >> "${OUT_DIR}/encode_comparison_${HISTONE}.csv"
}
 
# Initialize CSV if this is the first sample
if [ $SLURM_ARRAY_TASK_ID -eq 0 ]; then
    echo "Sample,Histone_Mark,Tool,Total_Peaks,ENCODE_Overlap,Recall_%,Precision_%" > "${OUT_DIR}/encode_comparison_H3K27ac.csv"
    echo "Sample,Histone_Mark,Tool,Total_Peaks,ENCODE_Overlap,Recall_%,Precision_%" > "${OUT_DIR}/encode_comparison_H3K27me3.csv"
fi
 
# Wait for CSV headers to be created
sleep 2
 
# Run comparisons for all 4 peak callers
calculate_overlap "MACS2" "${MACS2_DIR}/${SAMPLE}_dedup/${SAMPLE}_dedup_peaks.narrowPeak"
calculate_overlap "GoPeaks" "${GOPEAKS_DIR}/${SAMPLE}_dedup/${SAMPLE}_dedup_peaks.bed"
calculate_overlap "SEACR_Stringent" "${SEACR_DIR}/${SAMPLE}_fragments/${SAMPLE}_fragments.stringent.bed"
calculate_overlap "SEACR_Relaxed" "${SEACR_DIR}/${SAMPLE}_fragments/${SAMPLE}_fragments.relaxed.bed"

