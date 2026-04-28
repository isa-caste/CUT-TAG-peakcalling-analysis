#!/bin/bash
#SBATCH --mail-user=isacaste@iu.edu
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --job-name=bam-summary
#SBATCH -o bam-summary.out
#SBATCH -e bam-summary.err
#SBATCH --nodes=2
#SBATCH --mem=16g
#SBATCH -p gpu
#SBATCH --ntasks-per-node=2
#SBATCH --gpus-per-node=1
#SBATCH --time=1-23:59:00
#SBATCH -A r00750

set -e

QC_OUTDIR="/N/project/Krolab/isabella/data/bam-qc"
SUMMARY_FILE="$QC_OUTDIR/bam_qc_summary.csv"


# Process each flagstat file
for FLAGSTAT in "$QC_OUTDIR"/*_flagstat.txt; do
    if [ -f "$FLAGSTAT" ]; then
        BASENAME=$(basename "$FLAGSTAT" _flagstat.txt)
        
        # Extract stats from flagstat
        TOTAL_READS=$(head -1 "$FLAGSTAT" | awk '{print $1}')
        MAPPED_READS=$(grep "mapped (" "$FLAGSTAT" | head -1 | awk '{print $1}')
        
        if [ "$TOTAL_READS" -gt 0 ]; then
            PERCENT_MAPPED=$(echo "scale=2; ($MAPPED_READS / $TOTAL_READS) * 100" | bc)
        else
            PERCENT_MAPPED="0"
        fi
        
        # Get insert size stats
        INSERT_FILE="$QC_OUTDIR/${BASENAME}_insert_size_dist.txt"
        if [ -f "$INSERT_FILE" ]; then
            # Calculate mean insert size
            MEAN_INSERT=$(awk '{sum += $2 * $1; count += $1} END {if (count > 0) print int(sum/count)}' "$INSERT_FILE")
            MAX_INSERT=$(tail -1 "$INSERT_FILE" | awk '{print $2}')
        else
            MEAN_INSERT="NA"
            MAX_INSERT="NA"
        fi
        
        # Print CSV line (no units in the data, just numbers)
        echo "$BASENAME,$TOTAL_READS,$MAPPED_READS,$PERCENT_MAPPED,$MEAN_INSERT,$MAX_INSERT" >> "$SUMMARY_FILE"
    fi
done

echo "Summary report saved to: $SUMMARY_FILE"
