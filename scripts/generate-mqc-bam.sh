#!/bin/bash
#SBATCH --mail-user=isacaste@iu.edu
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH -o mqcbam.out
#SBATCH -e mqcbam.err
#SBATCH --job-name=mqcbam
#SBATCH -p gpu
#SBATCH --gpus-per-node=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=16G
#SBATCH --nodes=2
#SBATCH --time=1-23:59:00
#SBATCH -A r00750

set -e

BAM_DIR="/N/project/Krolab/isabella/data/bam-files"
QC_OUTDIR="/N/project/Krolab/isabella/data/bam-qc"
STATS_DIR="$QC_OUTDIR/samtools_stats"
MULTIQC_OUTDIR="$QC_OUTDIR/multiqc_report"

# Create output directory
mkdir -p "$MULTIQC_OUTDIR"

# Run MultiQC on samtools stats
echo "Running MultiQC on samtools stats..."
multiqc "$STATS_DIR" -o "$MULTIQC_OUTDIR" -n multiqc_report

echo "MultiQC Report Generated!"
