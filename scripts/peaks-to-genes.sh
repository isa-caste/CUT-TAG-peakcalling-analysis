#!/bin/bash
#SBATCH --job-name=peaks_to_genes
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/go-enrichment/logs/peaks_to_genes_%A_%a.out
#SBATCH --error=/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/go-enrichment/logs/peaks_to_genes_%A_%a.err
#SBATCH --array=0-1
#SBATCH --mail-user=isacaste@iu.edu
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH -A r00750
# Map peaks to nearest genes using bedtools
# Creates gene lists for GO enrichment

module load bedtools

# Paths
PEAK_DIR="/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-calling"
OUT_DIR="/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/go-enrichment"
GENES_FILE="/N/project/Krolab/isabella/data/annotations/hg38_genes.bed"

mkdir -p "${OUT_DIR}/logs"
mkdir -p "${OUT_DIR}/gene-lists"

# Peak caller options
TOOLS=("GoPeaks" "MACS2")
TOOL="${TOOLS[$SLURM_ARRAY_TASK_ID]}"

echo "Tool:       ${TOOL}"
echo "Array ID:   ${SLURM_ARRAY_TASK_ID}"
echo "============================="

# Check if gene annotation exists, if not download it
if [ ! -f "${GENES_FILE}" ]; then
    echo "Downloading hg38 gene annotations..."
    mkdir -p "$(dirname ${GENES_FILE})"
    
    # Download GENCODE genes
    wget -O - "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gtf.gz" | \
        gunzip | \
        awk '$3=="gene"' | \
        grep 'gene_type "protein_coding"' | \
        awk 'BEGIN{OFS="\t"}{
            match($0, /gene_id "([^"]+)"/, id);
            match($0, /gene_name "([^"]+)"/, name);
            print $1, $4-1, $5, name[1], ".", $7
        }' | \
        sort -k1,1 -k2,2n > "${GENES_FILE}"
    
    echo "Gene annotations downloaded: $(wc -l < ${GENES_FILE}) genes"
fi

# Representative samples (same as ChromHMM)
H3K27AC_SAMPLES=(SRR31972716 SRR31972718 SRR31972719 SRR31972720 SRR31972727)
H3K27ME3_SAMPLES=(SRR31972717 SRR31972726 SRR31972730 SRR31972735 SRR31972738 SRR31972739 SRR31972748)

ALL_SAMPLES=("${H3K27AC_SAMPLES[@]}" "${H3K27ME3_SAMPLES[@]}")

# Histone mark lookup
declare -A MARK_MAP
MARK_MAP["SRR31972716"]="H3K27ac"
MARK_MAP["SRR31972717"]="H3K27me3"
MARK_MAP["SRR31972718"]="H3K27ac"
MARK_MAP["SRR31972719"]="H3K27ac"
MARK_MAP["SRR31972720"]="H3K27ac"
MARK_MAP["SRR31972726"]="H3K27me3"
MARK_MAP["SRR31972727"]="H3K27ac"
MARK_MAP["SRR31972730"]="H3K27me3"
MARK_MAP["SRR31972735"]="H3K27me3"
MARK_MAP["SRR31972738"]="H3K27me3"
MARK_MAP["SRR31972739"]="H3K27me3"
MARK_MAP["SRR31972748"]="H3K27me3"

# Collect all genes across samples for each mark
H3K27AC_GENES="${OUT_DIR}/gene-lists/${TOOL}_H3K27ac_genes.txt"
H3K27ME3_GENES="${OUT_DIR}/gene-lists/${TOOL}_H3K27me3_genes.txt"

> "${H3K27AC_GENES}"
> "${H3K27ME3_GENES}"

# Process each sample
for SAMPLE in "${ALL_SAMPLES[@]}"; do
    
    MARK="${MARK_MAP[$SAMPLE]}"
    echo ""
    echo "Processing ${SAMPLE} (${MARK})..."
    
    # Get peak file
    if [ "${TOOL}" == "MACS2" ]; then
        PEAK_FILE="${PEAK_DIR}/macs2/${SAMPLE}_dedup/${SAMPLE}_dedup_peaks.narrowPeak"
    elif [ "${TOOL}" == "GoPeaks" ]; then
        PEAK_FILE="${PEAK_DIR}/gopeaks/${SAMPLE}_dedup/${SAMPLE}_dedup_peaks.bed"
    fi
    
    if [ ! -f "${PEAK_FILE}" ]; then
        echo "WARNING: Peak file not found: ${PEAK_FILE}"
        continue
    fi
    
    PEAK_COUNT=$(wc -l < "${PEAK_FILE}")
    echo "  Peaks: ${PEAK_COUNT}"
    
    # Find nearest gene for each peak
    NEAREST_FILE="${OUT_DIR}/gene-lists/${SAMPLE}_${TOOL}_nearest_genes.bed"
    
    bedtools closest \
        -a "${PEAK_FILE}" \
        -b "${GENES_FILE}" \
        -d \
        > "${NEAREST_FILE}"
    
    # Extract gene names (column 10 from bedtools closest output)
    cut -f10 "${NEAREST_FILE}" | sort -u >> "${OUT_DIR}/gene-lists/${TOOL}_${MARK}_genes_temp.txt"
    
    echo "Genes mapped"
    
done

# Combine and deduplicate genes for each mark
if [ -f "${OUT_DIR}/gene-lists/${TOOL}_H3K27ac_genes_temp.txt" ]; then
    sort -u "${OUT_DIR}/gene-lists/${TOOL}_H3K27ac_genes_temp.txt" > "${H3K27AC_GENES}"
    rm "${OUT_DIR}/gene-lists/${TOOL}_H3K27ac_genes_temp.txt"
    echo ""
    echo "H3K27ac gene list: $(wc -l < ${H3K27AC_GENES}) unique genes"
fi

if [ -f "${OUT_DIR}/gene-lists/${TOOL}_H3K27me3_genes_temp.txt" ]; then
    sort -u "${OUT_DIR}/gene-lists/${TOOL}_H3K27me3_genes_temp.txt" > "${H3K27ME3_GENES}"
    rm "${OUT_DIR}/gene-lists/${TOOL}_H3K27me3_genes_temp.txt"
    echo ""
    echo "H3K27me3 gene list: $(wc -l < ${H3K27ME3_GENES}) unique genes"
fi

echo " Gene mapping complete for ${TOOL}"
echo "============================="
