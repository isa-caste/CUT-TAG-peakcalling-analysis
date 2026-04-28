#!/bin/bash

# Check Duplication Rates Script
# Extracts duplication metrics from samtools stats files


set -e

STATS_DIR="/N/project/Krolab/isabella/data/bam-qc/samtools_stats"
OUTPUT_FILE="/N/project/Krolab/isabella/data/bam-qc/duplication_summary.txt"

echo "=========================================="
echo "Extracting Duplication Rates"
echo "=========================================="
echo ""

# Check if stats directory exists
if [ ! -d "$STATS_DIR" ]; then
    echo "ERROR: Stats directory not found: $STATS_DIR"
    echo "Please run bam_qc.sh first"
    exit 1
fi

# Count stats files
STATS_COUNT=$(ls "$STATS_DIR"/*_stats.txt 2>/dev/null | wc -l)
if [ $STATS_COUNT -eq 0 ]; then
    echo "ERROR: No stats files found in $STATS_DIR"
    echo "Please run bam_qc.sh first"
    exit 1
fi

echo "Found $STATS_COUNT stats files"
echo "Processing..."
echo ""

# Create summary file with header
echo "Sample Name | Mapped Reads | Duplicated Reads | Duplication Rate (%)" > "$OUTPUT_FILE"
echo "=================================================================" >> "$OUTPUT_FILE"

TOTAL_MAPPED=0
TOTAL_DUPLICATED=0
SAMPLE_COUNT=0
HIGH_DUP_SAMPLES=()
LOW_DUP_SAMPLES=()

# Process each stats file
for STATS_FILE in "$STATS_DIR"/*_stats.txt; do
    if [ -f "$STATS_FILE" ]; then
        BASENAME=$(basename "$STATS_FILE" _stats.txt)
        SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
        
        # Extract metrics from samtools stats
        # Look for lines like: SN	reads mapped and paired:	9400000
        MAPPED=$(grep "^SN.*reads mapped and paired" "$STATS_FILE" | awk '{print $NF}')
        DUPLICATED=$(grep "^SN.*reads duplicated" "$STATS_FILE" | awk '{print $NF}')
        
        # Handle cases where metrics might be missing
        if [ -z "$MAPPED" ]; then
            MAPPED=0
        fi
        if [ -z "$DUPLICATED" ]; then
            DUPLICATED=0
        fi
        
        # Calculate duplication rate
        if [ "$MAPPED" -gt 0 ]; then
            DUP_RATE=$(echo "scale=2; ($DUPLICATED / $MAPPED) * 100" | bc)
            TOTAL_MAPPED=$((TOTAL_MAPPED + MAPPED))
            TOTAL_DUPLICATED=$((TOTAL_DUPLICATED + DUPLICATED))
            
            # Store in arrays for categorization
            if (( $(echo "$DUP_RATE < 2" | bc -l) )); then
                LOW_DUP_SAMPLES+=("$BASENAME: $DUP_RATE%")
            elif (( $(echo "$DUP_RATE > 10" | bc -l) )); then
                HIGH_DUP_SAMPLES+=("$BASENAME: $DUP_RATE%")
            fi
            
            # Format output line
            printf "%-20s | %12s | %16s | %8s%%\n" \
                "$BASENAME" "$MAPPED" "$DUPLICATED" "$DUP_RATE" >> "$OUTPUT_FILE"
        fi
    fi
done

# Calculate overall statistics
if [ $SAMPLE_COUNT -gt 0 ] && [ $TOTAL_MAPPED -gt 0 ]; then
    OVERALL_DUP_RATE=$(echo "scale=2; ($TOTAL_DUPLICATED / $TOTAL_MAPPED) * 100" | bc)
else
    OVERALL_DUP_RATE=0
fi

# Add summary statistics
echo "" >> "$OUTPUT_FILE"
echo "SUMMARY STATISTICS" >> "$OUTPUT_FILE"
echo "=================================================================" >> "$OUTPUT_FILE"
echo "Total samples processed: $SAMPLE_COUNT" >> "$OUTPUT_FILE"
echo "Total mapped reads: $TOTAL_MAPPED" >> "$OUTPUT_FILE"
echo "Total duplicated reads: $TOTAL_DUPLICATED" >> "$OUTPUT_FILE"
echo "Overall duplication rate: $OVERALL_DUP_RATE%" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"

# Print to console
echo "=========================================="
echo "DUPLICATION RATE SUMMARY"
echo "=========================================="
echo ""
cat "$OUTPUT_FILE"
echo ""

# Categorize results
echo "=========================================="
echo "DATA QUALITY ASSESSMENT"
echo "=========================================="
echo ""

if (( $(echo "$OVERALL_DUP_RATE < 2" | bc -l) )); then
    echo "✓ EXCELLENT: Overall duplication rate is <2%"
    echo "  → Data quality is very high"
    echo "  → Duplicate removal is OPTIONAL (minimal benefit)"
    echo "  → But RECOMMENDED for consistency with paper"
elif (( $(echo "$OVERALL_DUP_RATE < 5" | bc -l) )); then
    echo "✓ GOOD: Overall duplication rate is 2-5%"
    echo "  → Data quality is good"
    echo "  → Duplicate removal is RECOMMENDED"
elif (( $(echo "$OVERALL_DUP_RATE < 10" | bc -l) )); then
    echo "⚠ ACCEPTABLE: Overall duplication rate is 5-10%"
    echo "  → Data quality is acceptable"
    echo "  → Duplicate removal is REQUIRED"
else
    echo "✗ HIGH: Overall duplication rate is >10%"
    echo "  → Data quality has issues"
    echo "  → Duplicate removal is REQUIRED"
    echo "  → Investigate potential problems"
fi

echo ""
echo "Expected for CUT&Tag: 1-5% duplicates"
echo ""

# Check for outliers
if [ ${#HIGH_DUP_SAMPLES[@]} -gt 0 ]; then
    echo "⚠ SAMPLES WITH HIGH DUPLICATION (>10%):"
    for SAMPLE in "${HIGH_DUP_SAMPLES[@]}"; do
        echo "  - $SAMPLE"
    done
    echo ""
fi

if [ ${#LOW_DUP_SAMPLES[@]} -gt 0 ]; then
    echo "✓ SAMPLES WITH EXCELLENT DATA (<2% duplication):"
    for SAMPLE in "${LOW_DUP_SAMPLES[@]}"; do
        echo "  - $SAMPLE"
    done
    echo ""
fi

echo "=========================================="
echo "RECOMMENDATION"
echo "=========================================="
echo ""
echo "For your project (CUT&Tag analysis with paper reproduction):"
echo "  → ALWAYS remove duplicates for consistency with methods"
echo "  → Use Picard (matches original paper)"
echo "  → Or use Sambamba (modern faster alternative)"
echo ""

# Print output location
echo "Full summary saved to: $OUTPUT_FILE"
echo ""
echo "Next step: Remove duplicates with Picard or Sambamba"
echo "=========================================="

