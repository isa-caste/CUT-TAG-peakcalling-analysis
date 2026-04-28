#!/bin/bash
# Peak Caller Comparison - Count peaks across all tools

# make paths
MACS2_DIR=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-calling/macs2
GOPEAKS_DIR=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-calling/gopeaks
SEACR_DIR=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-calling/seacr
METADATA=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/metadata.xlsx
OUT_DIR=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-comparison
mkdir -p "${OUT_DIR}"
OUTPUT="${OUT_DIR}/peak_counts_summary.csv"

# Create CSV header
echo "Sample,Histone_Mark,Antibody,DNA_Method,PCR_Cycles,MACS2,GoPeaks,SEACR_Stringent,SEACR_Relaxed" > "${OUTPUT}"
 
# Sample list (38 samples)
samples=(SRR31972716 SRR31972717 SRR31972718 SRR31972719 SRR31972720 \
         SRR31972721 SRR31972722 SRR31972723 SRR31972724 SRR31972725 \
         SRR31972726 SRR31972727 SRR31972728 SRR31972729 SRR31972730 \
         SRR31972731 SRR31972732 SRR31972733 SRR31972734 SRR31972735 \
         SRR31972736 SRR31972737 SRR31972738 SRR31972739 SRR31972740 \
         SRR31972741 SRR31972742 SRR31972743 SRR31972744 SRR31972745 \
         SRR31972746 SRR31972747 SRR31972748 SRR31972749 SRR31972750 \
         SRR31972751 SRR31972752 SRR31972753)
# Loop through each sample
for sample in "${samples[@]}"; do
    echo "Processing ${sample}..."
    
    # Count peaks for each tool
    macs2_count=$(wc -l < "${MACS2_DIR}/${sample}_dedup/${sample}_dedup_peaks.narrowPeak" 2>/dev/null || echo "0")
    gopeaks_count=$(wc -l < "${GOPEAKS_DIR}/${sample}_dedup/${sample}_dedup_peaks.bed" 2>/dev/null || echo "0")
    seacr_stringent=$(wc -l < "${SEACR_DIR}/${sample}_fragments/${sample}_fragments.stringent.bed" 2>/dev/null || echo "0")
    seacr_relaxed=$(wc -l < "${SEACR_DIR}/${sample}_fragments/${sample}_fragments.relaxed.bed" 2>/dev/null || echo "0")
    
    echo "${sample},NA,NA,NA,NA,${macs2_count},${gopeaks_count},${seacr_stringent},${seacr_relaxed}" >> "${OUTPUT}"
done
 
# For now, add placeholder metadata (will join with metadata file later)
    echo "${sample},NA,NA,NA,NA,${macs2_count},${gopeaks_count},${seacr_stringent},${seacr_relaxed}" >> "${OUTPUT}"
done

# Quick summary
echo "Total samples processed: ${#samples[@]}"
echo ""
echo "Average peak counts:"
awk -F',' 'NR>1 {m+=$6; g+=$7; ss+=$8; sr+=$9; n++} END {print "MACS2:          " m/n; print "GoPeaks:        " g/n; print "SEACR Stringent:" ss/n; print "SEACR Relaxed:  " sr/n}' "${OUTPUT}"

