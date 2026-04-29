# Benchmarking CUT&Tag Peak Callers and Genome Assemblies for Histone Modification Detection
## Overview
This project reproduces and extends the analysis from Abbasova et al. (2025), which demonstrated that CUT&Tag can recover up to ~54% of ENCODE ChIP-seq histone acetylation peaks. Three peak calling strategies were benchmarked MACS2, SEACR, and GoPeaks — using 38 K562 CUT&Tag samples profiling H3K27ac and H3K27me3 histone modifications, all aligned to hg38. <br>

Here it is evaluated if genome build choice (hg19 vs. hg38) and peak caller selection (MACS2, SEACR, GoPeaks) affect the biological conclusions drawn from CUT&Tag data. Peaks are validated against ENCODE ChIP-seq and evaluated through ENCODE recall/precision, ChromHMM chromatin state enrichment, Gene Ontology (GO) enrichment, and transcription factor motif analysis. <br>

Key findings: GoPeaks outperforms MACS2 for CUT&Tag data, achieving higher precision (80.6% vs 73.1%) and substantially better enrichment in functionally relevant chromatin states (up to 53× improvement).

## Data 
- **38 K562 CUT&Tag samples**: H3K27ac (n=31), H3K27me3 (n=7)
- Raw data available from NCBI SRA: accessions `SRR31972716`–`SRR31972753`
  - Download with: `bash scripts/download.sh` or manually via `prefetch` + `fasterq-dump`
- ENCODE ChIP-seq reference peaks used for validation:
  - H3K27ac: `ENCFF044JNJ` (lifted hg19 to hg38)
  - H3K27me3: `ENCFF000BXB` (lifted hg19 to hg38)
- ChromHMM reference: Roadmap Epigenomics K562 15-state model (E123, hg38)
> Raw FASTQ, BAM, and intermediate files are not stored in this repository due to size.
> All data can be reproduced from the SRA accessions above using the provided scripts.

## Pipeline Overview
```
Raw FASTQ
   │
   ├── Quality Control (FastQC, MultiQC)
   ├── Trimming (TrimGalore)
   ├── Alignment (Bowtie2, hg38)
   ├── Duplicate Removal (Sambamba)
   ├── BAM to BED conversion
   │
   ├── Peak Calling
   │     ├── GoPeaks (--broad)
   │     ├── MACS2 (--format BAMPE --qvalue 1e-5 --nolambda --nomodel)
   │     └── SEACR (stringent + relaxed modes)
   │
   └── Validation & Functional Analysis
         ├── ENCODE recall & precision (bedtools intersect)
         ├── ChromHMM enrichment (Roadmap E123)
         ├── GO enrichment (Enrichr API)
         └── Motif analysis (HOMER)
```
## Repository Structure 
```
CUT-TAG-peakcalling-analysis/
├── README.md
├── scripts/                        # All analysis scripts
└── results/
    ├── qc-reports/                 # FastQC and MultiQC HTML reports
    ├── delete-duplicates/          # Deduplication summaries
    ├── peak-calling/               # Peak count summaries per tool
    ├── encode-comparison/          # Recall/precision CSVs and plots
    ├── chromhmm-analysis/          # ChromHMM enrichment scores and heatmaps
    ├── go-enrichment-analysis/     # GO enrichment results
    └── motif-analysis/             # HOMER motif results per tool/mark
```
## Dependancies   
All tools were run on an HPC cluster using either module load or conda environments.

### Module-loaded Tools

| Tool | Version | Purpose |
|------|---------|---------|
| FastQC | 0.12.1 | Raw read QC |
| TrimGalore | 0.6.10 | Adapter trimming |
| Bowtie2 | 2.5.1 | Alignment to hg38 |
| samtools | 1.17 | BAM processing |
| bedtools | 2.31.0 | Peak overlap, BED operations |
| HOMER | 4.11.1 | Motif analysis |
| UCSC liftOver | 464 | hg19 → hg38 coordinate conversion |
| sra-toolkit | 3.0.5 | Raw data download from NCBI SRA |

### Conda Environments

| Environment | Tool | Version | Setup |
|-------------|------|---------|-------|
| `macs2_env` | MACS2 | 2.2.9.1 | `conda create -n macs2_env python=3.9 -y && pip install macs2` |
| `gopeaks_env` | GoPeaks | 1.0.0 | `conda create -n gopeaks_env -y && conda install -c bioconda gopeaks -y` |
| `seacr_env` | SEACR | 1.3 | `conda create -n seacr_env -y && conda install -c conda-forge bedtools r-base -y` |
| `seacr_env` | R | 4.5.3 | included above |
| base | MultiQC | 1.34.dev0 | `pip install multiqc` |

### SEACR Manual Installation
SEACR was not available as a conda package and had to be downloaded manually:
```bash
mkdir -p tools/SEACR
wget https://raw.githubusercontent.com/FredHutch/SEACR/master/SEACR_1.3.sh
wget https://raw.githubusercontent.com/FredHutch/SEACR/master/SEACR_1.3.R
chmod +x tools/SEACR/SEACR_1.3.sh
```
Update the `SEACR=` path in `scripts/seacr.sh` to point to your local copy.
## How to Reproduce 
1) Clone the repo
```
git clone https://github.com/<your-username>/CUT-TAG-peakcalling-analysis.git
cd CUT-TAG-peakcalling-analysis
```
2) Download the raw data
```
bash scripts/download.sh
```
3) Run QC and preprocessing
```
sbatch scripts/raw-fastqc-multiqc.sh
sbatch scripts/trimming.sh
sbatch scripts/align.sh
sbatch scripts/rm-dups.sh
```
4) Call Peaks with each tool
```
sbatch scripts/gopeaks-broad.sh
sbatch scripts/macs2.sh
sbatch scripts/seacr.sh
```
5) Run Validation & analyses
```
sbatch scripts/encode-comparison.sh
sbatch scripts/chromhmm-enrich.sh
python scripts/go-enrich.py
sbatch scripts/motif.sh
```
> NOTE <br>
> Scripts assume an SLURM HPC environment. Adjust resource parameters as needed for your system.
## References
1) Abbasova, L., Urbanaviciute, P., Hu, D. et al. CUT&Tag recovers up to half of ENCODE ChIP-seq histone acetylation peaks. Nat Commun 16, 2993 (2025). https://doi.org/10.1038/s41467-025-58137-2 <br>
2) Yashar, W. M., Kong, G., VanCampen, J., et al. GoPeaks: histone modification peak calling for CUT&Tag. Genome Biology 23, 144 (2022). https://doi.org/10.1186/s13059-022-02707-w
