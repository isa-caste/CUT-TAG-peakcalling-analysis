#!/bin/bash
# Download K562 Chromatin State Annotations (ChromHMM 15-state model)
# From Roadmap Epigenomics Consortium
# Cell line: E123 (K562)

# --- Directories ---
CHROMHMM_DIR="/N/project/Krolab/isabella/data/chromhmm"
mkdir -p "${CHROMHMM_DIR}"

cd "${CHROMHMM_DIR}"

echo "========================================="
echo "Downloading K562 Chromatin States (hg38)"
echo "========================================="

# K562 = E123 in Roadmap Epigenomics
CELL_TYPE="E123"

# Download 15-state ChromHMM model (hg38 liftOver)
# This file has already been lifted from hg19 to hg38
if [ ! -f "${CELL_TYPE}_15_coreMarks_hg38lift_mnemonics.bed.gz" ]; then
    echo "Downloading ${CELL_TYPE} chromatin states..."
    wget --no-check-certificate \
        https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/${CELL_TYPE}_15_coreMarks_hg38lift_mnemonics.bed.gz
    
    if [ $? -ne 0 ]; then
        echo "ERROR: Download failed. Trying alternative URL..."
        # Alternative: download hg19 version and we'll liftOver ourselves
        wget http://www.broadinstitute.org/~jernst/CHMM_15state_hg19/${CELL_TYPE}_15_coreMarks_mnemonics.bed.gz
    fi
else
    echo "ChromHMM file already exists, skipping download."
fi

# Uncompress
if [ -f "${CELL_TYPE}_15_coreMarks_hg38lift_mnemonics.bed.gz" ]; then
    echo "Uncompressing..."
    gunzip -f ${CELL_TYPE}_15_coreMarks_hg38lift_mnemonics.bed.gz
    CHROMHMM_BED="${CELL_TYPE}_15_coreMarks_hg38lift_mnemonics.bed"
elif [ -f "${CELL_TYPE}_15_coreMarks_mnemonics.bed.gz" ]; then
    echo "Uncompressing hg19 version..."
    gunzip -f ${CELL_TYPE}_15_coreMarks_mnemonics.bed.gz
    CHROMHMM_BED="${CELL_TYPE}_15_coreMarks_mnemonics.bed"
    echo "WARNING: This is hg19. You may need to liftOver to hg38."
else
    echo "ERROR: Download failed"
    exit 1
fi

echo ""
echo "========================================="
echo "File Information"
echo "========================================="

echo "File: ${CHROMHMM_BED}"
echo "Lines: $(wc -l < ${CHROMHMM_BED})"
echo "File size: $(du -h ${CHROMHMM_BED} | cut -f1)"

echo ""
echo "Chromatin state distribution:"
cut -f4 ${CHROMHMM_BED} | sort | uniq -c | sort -rn

echo ""
echo "Sample regions:"
head -10 ${CHROMHMM_BED}

echo ""
echo "========================================="
echo "Chromatin State Legend (15-state model)"
echo "========================================="
cat << 'EOF'

State  Name         Color        Biology
-----  -----------  -----------  ---------------------------
1      TssA         Red          Active TSS
2      TssFlnk      OrangeRed    Flanking TSS
3      TssFlnkU     OrangeRed    Flanking TSS Upstream
4      TssFlnkD     OrangeRed    Flanking TSS Downstream
5      Tx           Green        Strong transcription
6      TxWk         DarkGreen    Weak transcription
7      EnhG1        YellowGreen  Genic enhancer 1
8      EnhG2        YellowGreen  Genic enhancer 2
9      EnhA1        Orange       Active enhancer 1
10     EnhA2        Orange       Active enhancer 2
11     EnhWk        Yellow       Weak enhancer
12     ZNF/Rpts     Teal         ZNF genes & repeats
13     Het          Gray         Heterochromatin
14     TssBiv       DarkOrange   Bivalent/poised TSS
15     EnhBiv       LightOrange  Bivalent enhancer
16     ReprPC       Purple       Repressed PolyComb
17     ReprPCWk     LightPurple  Weak Repressed PolyComb
18     Quies        White        Quiescent/Low signal

EOF

echo "Download complete!"
echo "ChromHMM annotation file: ${CHROMHMM_DIR}/${CHROMHMM_BED}"

