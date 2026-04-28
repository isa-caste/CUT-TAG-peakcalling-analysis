#!/bin/bash
set -euo pipefail

# Validate deduplicated BAM files from Picard MarkDuplicates

module load conda
conda activate chip_env

DEDUP_DIR="/N/project/Krolab/isabella/data/dedup-bam-files"
OUTPUT_DIR="/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/dedup-validation"

mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "Validating Deduplicated BAM Files"
echo "=========================================="
echo "Input directory: $DEDUP_DIR"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Initialize summary CSV
SUMMARY_FILE="$OUTPUT_DIR/dedup_validation_summary.csv"
cat > "$SUMMARY_FILE" << 'EOF'
sample_name,read_pairs_examined,duplicates_marked,percent_dup,estimated_lib_size,bam_size_mb,status
EOF

# Process each sample
for metrics_file in "$DEDUP_DIR"/*_dedup_metrics.txt; do
    if [[ ! -f "$metrics_file" ]]; then
        continue
    fi

    sample_name=$(basename "$metrics_file" _dedup_metrics.txt)
    bam_file="${DEDUP_DIR}/${sample_name}_dedup.bam"
    
    echo "Sample: $sample_name"

    # Check if BAM exists
    if [[ ! -f "$bam_file" ]]; then
        echo "BAM file not found"
        echo "$sample_name,NA,NA,NA,NA,NA,MISSING_BAM" >> "$SUMMARY_FILE"
        continue
    fi

    # Get BAM file size in MB
    bam_size=$(du -m "$bam_file" | awk '{print $1}')

    # Parse Picard metrics file
    # Format: After header lines, look for the data line with LIBRARY and PERCENT_DUPLICATION
    read_pairs_examined=$(grep -E "^lib[0-9]|^LIBRARY" "$metrics_file" | grep -v "^LIBRARY" | awk '{print $3}')
    duplicates_marked=$(grep -E "^lib[0-9]|^LIBRARY" "$metrics_file" | grep -v "^LIBRARY" | awk '{print $7}')
    percent_dup=$(grep -E "^lib[0-9]|^LIBRARY" "$metrics_file" | grep -v "^LIBRARY" | awk '{print $9}')
    estimated_lib=$(grep -E "^lib[0-9]|^LIBRARY" "$metrics_file" | grep -v "^LIBRARY" | awk '{print $10}')

    if [[ -n "$read_pairs_examined" ]]; then
        # Convert decimal to percentage
        percent_dup_pct=$(echo "scale=2; $percent_dup * 100" | bc)
        echo "Read pairs examined: $read_pairs_examined"
        echo "Duplicates marked: $duplicates_marked"
        echo "Duplication rate: $percent_dup_pct%"
        echo "Estimated library size: $estimated_lib"
        echo "BAM size: ${bam_size}M"
        echo "PASS"
        status="PASS"
    else
        echo "Could not parse metrics file"
        read_pairs_examined="NA"
        duplicates_marked="NA"
        percent_dup_pct="NA"
        estimated_lib="NA"
        status="PARSE_ERROR"
    fi

    echo ""
    
    # Append to summary
    echo "$sample_name,$read_pairs_examined,$duplicates_marked,$percent_dup_pct,$estimated_lib,$bam_size,$status" >> "$SUMMARY_FILE"
done
