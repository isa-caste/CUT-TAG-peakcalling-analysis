#!/bin/bash
#SBATCH --mail-user=isacaste@iu.edu 
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH -o bamqc.out
#SBATCH -e bamqc.err
#SBATCH --job-name=bamqc
#SBATCH -p gpu
#SBATCH --nodes=2
#SBATCH --mem=16g
#SBATCH --ntasks-per-node=2
#SBATCH --gpus-per-node=1
#SBATCH --time=1-23:59:00
#SBATCH -A r00750

set -e

# load modules
module load samtools

# Define directories
BAM_DIR="/N/project/Krolab/isabella/data/bam-files"
QC_OUTDIR="/N/project/Krolab/isabella/data/bam-qc"
STATS_DIR="$QC_OUTDIR/samtools_stats"

# Create output directories
mkdir -p "$QC_OUTDIR"
mkdir -p "$STATS_DIR"

# Counter for tracking progress
TOTAL=$(ls "$BAM_DIR"/*_sorted.bam 2>/dev/null | wc -l)
COUNT=0

# Process each BAM file
for BAM in "$BAM_DIR"/*_sorted.bam; do
    COUNT=$((COUNT + 1))
    BASENAME=$(basename "$BAM" _sorted.bam)
    
    echo "[$COUNT/$TOTAL] Processing $BASENAME..."
    
    # 1. Flagstat - Basic alignment statistics
    echo "  - Running flagstat..."
    samtools flagstat "$BAM" > "$QC_OUTDIR/${BASENAME}_flagstat.txt"
    
    # 2. Samtools stats - Detailed statistics
    echo "  - Running samtools stats..."
    samtools stats "$BAM" > "$STATS_DIR/${BASENAME}_stats.txt"
    
    # 3. Insert size distribution (for paired-end reads)
    echo "  - Calculating insert size distribution..."
    samtools view -f 0x2 "$BAM" | \
      awk '{if ($9 > 0) print $9}' | \
      sort -n | uniq -c > "$QC_OUTDIR/${BASENAME}_insert_size_dist.txt"
    
    # 4. Count total reads
    echo "  - Counting reads..."
    TOTAL_READS=$(samtools view -c "$BAM")
    MAPPED_READS=$(samtools view -c -F 4 "$BAM")
    PERCENT_MAPPED=$(echo "scale=2; ($MAPPED_READS / $TOTAL_READS) * 100" | bc)
    
    echo "    Total reads: $TOTAL_READS"
    echo "    Mapped reads: $MAPPED_READS ($PERCENT_MAPPED%)"
    
    # 5. Fragment length statistics
    echo "  - Calculating fragment length statistics..."
    MEAN_INSERT=$(samtools view -f 0x2 "$BAM" | awk '{sum+=$9; count++} END {if (count>0) print int(sum/count)}')
    echo "    Mean insert size: $MEAN_INSERT bp"
    
done

echo "BAM QC Complete!"
echo ""
echo "Output files generated:"
echo "  - Flagstat reports: $QC_OUTDIR/*_flagstat.txt"
echo "  - Detailed stats: $STATS_DIR/*_stats.txt"
echo "  - Insert size distribution: $QC_OUTDIR/*_insert_size_dist.txt"
