#!/bin/bash
# Download ENCODE reference files and liftOver from hg19 to hg38
# Fixed version: uses ucsctools module and converts to BED3 format before liftOver

# --- Directories ---
ENCODE_DIR="/N/project/Krolab/isabella/data/encode-files"
LIFTOVER_DIR="${ENCODE_DIR}/hg38"
LOGS_DIR="${ENCODE_DIR}/logs"

mkdir -p "${ENCODE_DIR}" "${LIFTOVER_DIR}" "${LOGS_DIR}"

echo "========================================="
echo "Step 1: Download ENCODE files"
echo "========================================="

cd "${ENCODE_DIR}"

# Blacklist regions (bigBed)
if [ ! -f "ENCFF000KJP.bigBed" ]; then
    echo "Downloading ENCFF000KJP (blacklist)..."
    wget -q --show-progress \
        https://www.encodeproject.org/files/ENCFF000KJP/@@download/ENCFF000KJP.bigBed
else
    echo "ENCFF000KJP.bigBed already exists, skipping."
fi

# H3K27me3 ChIP-seq peaks (bigBed)
if [ ! -f "ENCFF000BXB.bigBed" ]; then
    echo "Downloading ENCFF000BXB (H3K27me3 peaks)..."
    wget -q --show-progress \
        https://www.encodeproject.org/files/ENCFF000BXB/@@download/ENCFF000BXB.bigBed
else
    echo "ENCFF000BXB.bigBed already exists, skipping."
fi

# H3K27ac ChIP-seq peaks (bed.gz)
if [ ! -f "ENCFF044JNJ.bed" ]; then
    echo "Downloading ENCFF044JNJ (H3K27ac peaks)..."
    wget -q --show-progress \
        https://www.encodeproject.org/files/ENCFF044JNJ/@@download/ENCFF044JNJ.bed.gz
    gunzip ENCFF044JNJ.bed.gz
else
    echo "ENCFF044JNJ.bed already exists, skipping."
fi

echo ""
echo "========================================="
echo "Step 2: Convert bigBed → BED"
echo "========================================="

# Load UCSC tools module (correct module name for this HPC)
module load ucsctools/464

# Check bigBedToBed is available
if ! command -v bigBedToBed &> /dev/null; then
    echo "ERROR: bigBedToBed not found after loading ucsctools/464"
    exit 1
fi

echo "Converting ENCFF000KJP.bigBed → ENCFF000KJP.bed..."
bigBedToBed ENCFF000KJP.bigBed ENCFF000KJP.bed

echo "Converting ENCFF000BXB.bigBed → ENCFF000BXB.bed..."
bigBedToBed ENCFF000BXB.bigBed ENCFF000BXB.bed

echo ""
echo "========================================="
echo "Step 3: Download hg19 → hg38 liftOver chain file"
echo "========================================="

CHAIN_FILE="${ENCODE_DIR}/hg19ToHg38.over.chain.gz"

if [ ! -f "${CHAIN_FILE}" ]; then
    echo "Downloading liftOver chain file..."
    wget -q --show-progress \
        https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
else
    echo "Chain file already exists, skipping."
fi

echo ""
echo "========================================="
echo "Step 4: Convert to BED3 format (required for liftOver)"
echo "========================================="

# liftOver cannot handle narrowPeak/broadPeak format (10 columns)
# Extract only chr, start, end (first 3 columns)
echo "Converting ENCFF000KJP to BED3..."
cut -f1-3 ENCFF000KJP.bed > ENCFF000KJP_bed3.bed

echo "Converting ENCFF000BXB to BED3..."
cut -f1-3 ENCFF000BXB.bed > ENCFF000BXB_bed3.bed

echo "Converting ENCFF044JNJ to BED3..."
cut -f1-3 ENCFF044JNJ.bed > ENCFF044JNJ_bed3.bed

echo ""
echo "========================================="
echo "Step 5: LiftOver all files to hg38"
echo "========================================="

# ucsctools module already loaded above, includes liftOver

if ! command -v liftOver &> /dev/null; then
    echo "ERROR: liftOver not found after loading ucsctools/464"
    exit 1
fi

# liftOver blacklist
echo "LiftOver: ENCFF000KJP (blacklist) hg19 → hg38..."
liftOver \
    ENCFF000KJP_bed3.bed \
    "${CHAIN_FILE}" \
    "${LIFTOVER_DIR}/ENCFF000KJP_hg38.bed" \
    "${LIFTOVER_DIR}/ENCFF000KJP_unmapped.bed"

# liftOver H3K27me3 peaks
echo "LiftOver: ENCFF000BXB (H3K27me3) hg19 → hg38..."
liftOver \
    ENCFF000BXB_bed3.bed \
    "${CHAIN_FILE}" \
    "${LIFTOVER_DIR}/ENCFF000BXB_hg38.bed" \
    "${LIFTOVER_DIR}/ENCFF000BXB_unmapped.bed"

# liftOver H3K27ac peaks
echo "LiftOver: ENCFF044JNJ (H3K27ac) hg19 → hg38..."
liftOver \
    ENCFF044JNJ_bed3.bed \
    "${CHAIN_FILE}" \
    "${LIFTOVER_DIR}/ENCFF044JNJ_hg38.bed" \
    "${LIFTOVER_DIR}/ENCFF044JNJ_unmapped.bed"

echo ""
echo "========================================="
echo "Step 6: Summary"
echo "========================================="

echo "File                    Mapped    Unmapped"
echo "-------------------------------------------"

for f in "${LIFTOVER_DIR}"/*_hg38.bed; do
    name=$(basename "$f" _hg38.bed)
    mapped=$(wc -l < "$f")
    unmapped_file="${LIFTOVER_DIR}/${name}_unmapped.bed"
    unmapped=$(wc -l < "$unmapped_file" 2>/dev/null || echo "0")
    printf "%-20s %8d %10d\n" "$name" "$mapped" "$unmapped"
done

echo ""
echo "✓ Done! hg38 files are in: ${LIFTOVER_DIR}"
echo ""
echo "ENCODE reference peaks (hg38):"
echo "  H3K27ac:      ${LIFTOVER_DIR}/ENCFF044JNJ_hg38.bed ($(wc -l < ${LIFTOVER_DIR}/ENCFF044JNJ_hg38.bed) peaks)"
echo "  H3K27me3:     ${LIFTOVER_DIR}/ENCFF000BXB_hg38.bed ($(wc -l < ${LIFTOVER_DIR}/ENCFF000BXB_hg38.bed) peaks)"
echo "  Blacklist:    ${LIFTOVER_DIR}/ENCFF000KJP_hg38.bed ($(wc -l < ${LIFTOVER_DIR}/ENCFF000KJP_hg38.bed) regions)"
