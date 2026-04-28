#!/bin/bash
#SBATCH --account=r00750
set -euo pipefail

#SBATCH --job-name=rem-dups
#SBATCH --mail-user=isacaste@iu.edu
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH -o rem-dups.out
#SBATCH -e rem-dups.err
#SBATCH -p h2
#SBATCH --ntasks=8
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH -A r00750

set +u
module load conda
conda activate chip_env
set -u

# Verify Picard is available
if ! command -v picard &> /dev/null; then
    echo "ERROR: picard not found in chip_env. Installing..."
    conda install -c bioconda picard -y
fi

# Verify samtools is available
if ! command -v samtools &> /dev/null; then
    echo "ERROR: samtools not found. Please ensure it's in chip_env"
    exit 1
fi

BAM_DIR=/N/project/Krolab/isabella/data/bam-files
OUT_DIR=/N/project/Krolab/isabella/data/dedup-bam-files
LOG_DIR=/N/project/Krolab/isabella/data/dedup-bam-files/logs

# Create output directories
mkdir -p "$OUT_DIR" "$LOG_DIR"

echo "=========================================="
echo "Picard MarkDuplicates - Batch Processing"
echo "=========================================="
echo "Input BAM directory: $BAM_DIR"
echo "Output directory: $OUT_DIR"
echo "Timestamp: $(date)"
echo ""

# Count input BAMs
BAM_COUNT=$(ls -1 "$BAM_DIR"/*_sorted.bam 2>/dev/null | wc -l)
echo "Found $BAM_COUNT sorted BAM files to process"
echo ""

PROCESSED=0
FAILED=0

for BAM in "$BAM_DIR"/*_sorted.bam; do
    if [[ ! -f "$BAM" ]]; then
        continue
    fi

    SAMPLE=$(basename "$BAM" _sorted.bam)
    SAMPLE_LOG="$LOG_DIR/${SAMPLE}.log"

    echo "======== Processing: $SAMPLE ========"
    {
        echo "Start time: $(date)"
        echo "BAM input: $BAM"
        
        # Step 1: Verify input BAM exists and is readable
        if ! samtools quickcheck "$BAM" 2>&1; then
            echo "ERROR: Input BAM file is corrupted or unreadable"
            exit 1
        fi
        echo " Input BAM validated"

        # Step 2: Add or replace read groups
        echo "Step 1/3: Adding read groups..."
        RG_BAM="$OUT_DIR/${SAMPLE}_rg.bam"
        picard AddOrReplaceReadGroups \
            I="$BAM" \
            O="$RG_BAM" \
            RGID="$SAMPLE" \
            RGLB="lib1" \
            RGPL="illumina" \
            RGPU="unit1" \
            RGSM="$SAMPLE" \
            VALIDATION_STRINGENCY=LENIENT \
            2>&1 | tail -20
        
        if [[ ! -f "$RG_BAM" ]]; then
            echo "ERROR: Read group BAM not created"
            exit 1
        fi
        echo " Read groups added: $RG_BAM"

        # Step 3: Mark and remove duplicates
        echo "Step 2/3: Marking and removing duplicates..."
        DEDUP_BAM="$OUT_DIR/${SAMPLE}_dedup.bam"
        METRICS_FILE="$OUT_DIR/${SAMPLE}_dedup_metrics.txt"
        
        picard MarkDuplicates \
            I="$RG_BAM" \
            O="$DEDUP_BAM" \
            M="$METRICS_FILE" \
            REMOVE_DUPLICATES=true \
            ASSUME_SORTED=true \
            VALIDATION_STRINGENCY=LENIENT \
            2>&1 | tail -20
        
        if [[ ! -f "$DEDUP_BAM" ]]; then
            echo "ERROR: Deduplicated BAM not created"
            exit 1
        fi
        echo " Duplicates removed: $DEDUP_BAM"

        # Step 4: Index the deduplicated BAM
        echo "Step 3/3: Creating BAM index..."
        samtools index -@ 8 "$DEDUP_BAM"
        
        if [[ ! -f "${DEDUP_BAM}.bai" ]]; then
            echo "ERROR: Index not created"
            exit 1
        fi
        echo " Index created: ${DEDUP_BAM}.bai"

        # Step 5: Cleanup intermediate file
        echo "Cleaning up intermediate files..."
        rm -f "$RG_BAM"
        echo "Intermediate files removed"

        # Step 6: Verify outputs
        echo "Verifying final outputs..."
        if samtools quickcheck "$DEDUP_BAM" 2>&1; then
            echo "Final BAM validation passed"
        else
            echo "WARNING: Final BAM validation failed"
            exit 1
        fi

        # Print duplication statistics
        echo ""
        echo "Duplication Statistics:"
        if [[ -f "$METRICS_FILE" ]]; then
            # Skip header lines and print data
            tail -n +8 "$METRICS_FILE" | head -1
        fi

        echo " Sample $SAMPLE completed successfully"
        echo "End time: $(date)"

    } 2>&1 | tee -a "$SAMPLE_LOG"

    if [[ ${PIPESTATUS[0]} -eq 0 ]]; then
        ((PROCESSED++))
        echo " $SAMPLE: SUCCESS"
    else
        ((FAILED++))
        echo " $SAMPLE: FAILED (see $SAMPLE_LOG)"
    fi
    echo ""
done

echo "=========================================="
echo "SUMMARY"
echo "=========================================="
echo "Total processed: $PROCESSED / $BAM_COUNT"
echo "Failed: $FAILED"
echo "Timestamp: $(date)"
echo ""

if [[ $FAILED -eq 0 ]]; then
    echo "All samples processed successfully!"
    exit 0
else
    echo " $FAILED sample(s) failed. Check logs in $LOG_DIR"
    exit 1
fi
