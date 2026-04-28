#!/bin/bash
#SBATCH --job-name=bam-to-bed-parallel
#SBATCH --mail-user=isacaste@iu.edu
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH -o bam-to-bed-parallel.out
#SBATCH -e bam-to-bed-parallel.err
#SBATCH -p h2
#SBATCH --ntasks=8
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --time=23:59:00
#SBATCH -A r00750

set -euo pipefail

set +u
module load conda
conda activate chip_env
set -u

# Verify bedtools is available
if ! command -v bedtools &> /dev/null; then
    echo "Installing bedtools..."
    conda install -c bioconda bedtools -y
fi

BAM_DIR=/N/project/Krolab/isabella/data/dedup-bam-files
BED_DIR=/N/project/Krolab/isabella/data/bed-files

mkdir -p "$BED_DIR"

echo "=========================================="
echo "Converting BAM to BED (parallel)"
echo "=========================================="
echo "Input BAM directory: $BAM_DIR"
echo "Output BED directory: $BED_DIR"
echo "Timestamp: $(date)"
echo ""

# Count input BAMs
BAM_COUNT=$(ls -1 "$BAM_DIR"/*_dedup.bam 2>/dev/null | wc -l)
echo "Found $BAM_COUNT deduplicated BAM files to process"
echo ""

# Create a function to process each BAM
process_bam() {
    BAM="$1"
    BED_DIR="$2"
    
    SAMPLE=$(basename "$BAM" _dedup.bam)
    BED_FILE="$BED_DIR/${SAMPLE}_fragments.bed"
    
    if bedtools bamtobed -i "$BAM" -bedpe 2>/dev/null | \
        awk 'BEGIN {OFS="\t"} {print $1, $2, $6}' | \
        sort -k1,1 -k2,2n > "$BED_FILE"; then
        
        BED_LINES=$(wc -l < "$BED_FILE")
        echo "$SAMPLE: $BED_LINES fragments ✓"
    else
        echo "$SAMPLE: FAILED ❌"
        return 1
    fi
}

export -f process_bam
export BED_DIR

# Process all BAMs in parallel (8 at a time, matching number of CPUs)
echo "Processing in parallel..."
ls "$BAM_DIR"/*_dedup.bam | \
    parallel -j 8 process_bam {} "$BED_DIR"

echo ""
echo "=========================================="
echo "SUMMARY"
echo "=========================================="
echo "Timestamp: $(date)"

# Count results
BED_COUNT=$(ls -1 "$BED_DIR"/*.bed 2>/dev/null | wc -l)
echo "Successfully created: $BED_COUNT / $BAM_COUNT BED files"

if [[ $BED_COUNT -eq $BAM_COUNT ]]; then
    echo "All samples processed!"
else
    echo "Only $BED_COUNT/$BAM_COUNT completed. Check error log."
fi

echo ""
echo "Output location: $BED_DIR"
echo "Ready for peak calling with MACS2, SEACR, or GoPeaks"
