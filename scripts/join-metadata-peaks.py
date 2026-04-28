#!/usr/bin/env python3
"""
Join peak counts with sample metadata
Adds histone mark, antibody, DNA extraction method, PCR cycles to peak count table
"""
import pandas as pd
import sys
 
# Paths
METADATA = '/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/metadata.xlsx'
PEAK_COUNTS = '/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-comparison/peak_counts_summary.csv'
OUTPUT = '/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-comparison/peak_counts_with_metadata.csv'
# Read metadata
print(f"\nReading metadata from: {METADATA}")
metadata = pd.read_excel(METADATA)
 
# Read peak counts
print(f"Reading peak counts from: {PEAK_COUNTS}")
peaks = pd.read_csv(PEAK_COUNTS)
 
# Clean up metadata - keep only relevant columns
metadata_clean = metadata[['Run', 'histone_mark', 'antibody', 'DNA_extraction_method', 'PCR_cycles']].copy()
# Simplify antibody names
metadata_clean['antibody_simple'] = (
    metadata_clean['antibody']
    .str.replace('Diagenode C15410196', 'Diagenode')
    .str.replace('Active Motif 39133', 'ActiveMotif')
    .str.replace('Abcam-ab4729', 'Abcam-4729')
    .str.replace('Abcam-ab177178', 'Abcam-177178')
    .str.replace('CST 9733', 'CST')
)
# Merge
print("\nMerging datasets...")
merged = peaks.merge(
    metadata_clean[['Run', 'histone_mark', 'antibody_simple', 'DNA_extraction_method', 'PCR_cycles']], 
    left_on='Sample', 
    right_on='Run', 
    how='left'
)
 
# Drop old placeholder columns and Run duplicate
merged = merged.drop(columns=['Histone_Mark', 'Antibody', 'DNA_Method', 'PCR_Cycles', 'Run'])
 
# Rename columns
merged = merged.rename(columns={
    'histone_mark': 'Histone_Mark',
    'antibody_simple': 'Antibody',
    'DNA_extraction_method': 'DNA_Method',
    'PCR_cycles': 'PCR_Cycles'
})
 
# Reorder columns
merged = merged[['Sample', 'Histone_Mark', 'Antibody', 'DNA_Method', 'PCR_Cycles', 
                 'MACS2', 'GoPeaks', 'SEACR_Stringent', 'SEACR_Relaxed']]
# Save
print(f"\nSaving to: {OUTPUT}")
merged.to_csv(OUTPUT, index=False)

print(f"\n Merged file created with {len(merged)} samples")
print("\nSample preview:")
print(merged.head(10).to_string())

print("Summary by Histone Mark:")
print("=" * 80)
summary = merged.groupby('Histone_Mark')[['MACS2', 'GoPeaks', 'SEACR_Stringent', 'SEACR_Relaxed']].mean()
print(summary)

print("Summary by Antibody (H3K27ac only):")
print("=" * 80)
h3k27ac = merged[merged['Histone_Mark'] == 'H3K27ac']
summary_ab = h3k27ac.groupby('Antibody')[['MACS2', 'GoPeaks', 'SEACR_Stringent', 'SEACR_Relaxed']].mean()
print(summary_ab)

