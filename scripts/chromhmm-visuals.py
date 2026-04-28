#!/usr/bin/env python3
"""
Visualize ChromHMM State Enrichment
Creates heatmap showing enrichment of peaks in different chromatin states
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

# Paths
GOPEAKS_FILE = '/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/chromhmm/GoPeaks_chromhmm_enrichment.csv'
MACS2_FILE = '/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/chromhmm/MACS2_chromhmm_enrichment.csv'
OUT_DIR = '/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/chromhmm/figures'

os.makedirs(OUT_DIR, exist_ok=True)

# Set plotting style
sns.set_style("white")
plt.rcParams['font.size'] = 10

print("=" * 80)
print("ChromHMM Enrichment Visualization")
print("=" * 80)

# Read data
print("\nReading enrichment data...")
gopeaks = pd.read_csv(GOPEAKS_FILE)
macs2 = pd.read_csv(MACS2_FILE)

# Combine
all_data = pd.concat([gopeaks, macs2], ignore_index=True)

print(f"Total records: {len(all_data)}")
print(f"  GoPeaks: {len(gopeaks)}")
print(f"  MACS2: {len(macs2)}")

# Simplify state names (remove numbers)
all_data['State_Simple'] = all_data['State'].str.replace(r'^\d+_', '', regex=True)

# Group states into categories for cleaner visualization
STATE_GROUPS = {
    'Active Promoter': ['TssA', 'TssFlnk', 'TssFlnkU', 'TssFlnkD'],
    'Transcription': ['Tx', 'TxWk', 'TxFlnk'],
    'Active Enhancer': ['Enh', 'EnhG'],
    'Bivalent': ['TssBiv', 'EnhBiv', 'BivFlnk'],
    'Repressed': ['ReprPC', 'ReprPCWk', 'Het'],
    'Low Signal': ['Quies', 'ZNF/Rpts']
}

# Add state category
def categorize_state(state):
    for category, states in STATE_GROUPS.items():
        if any(s in state for s in states):
            return category
    return 'Other'

all_data['State_Category'] = all_data['State_Simple'].apply(categorize_state)

# Average enrichment by tool, histone mark, and state category
print("\nCalculating average enrichment scores...")
enrichment_summary = all_data.groupby(['Tool', 'Histone_Mark', 'State_Category'])['Enrichment_Score'].mean().reset_index()

print("\n" + "=" * 80)
print("ENRICHMENT SUMMARY")
print("=" * 80)
print(enrichment_summary.pivot_table(index='State_Category', 
                                      columns=['Histone_Mark', 'Tool'], 
                                      values='Enrichment_Score'))

# ===========================
# Figure 1: Heatmap of enrichment by state category
# ===========================
print("\nGenerating heatmap...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for idx, mark in enumerate(['H3K27ac', 'H3K27me3']):
    subset = enrichment_summary[enrichment_summary['Histone_Mark'] == mark]
    
    # Pivot for heatmap
    pivot = subset.pivot(index='State_Category', columns='Tool', values='Enrichment_Score')
    
    # Reorder rows by biological relevance
    row_order = ['Active Promoter', 'Active Enhancer', 
                 'Transcription', 'Bivalent', 'Repressed', 'Low Signal']
    pivot = pivot.reindex([r for r in row_order if r in pivot.index])
    
    # Create heatmap
    sns.heatmap(pivot, 
                annot=True, 
                fmt='.2f', 
                cmap='RdYlGn', 
                center=1.0,  # Enrichment of 1 = expected
                vmin=0, 
                vmax=3,
                cbar_kws={'label': 'Enrichment Score'},
                ax=axes[idx])
    
    axes[idx].set_title(f'{mark} - Chromatin State Enrichment', fontsize=14, fontweight='bold')
    axes[idx].set_xlabel('Peak Caller', fontsize=12)
    axes[idx].set_ylabel('Chromatin State Category', fontsize=12)

plt.tight_layout()
plt.savefig(f'{OUT_DIR}/chromhmm_enrichment_heatmap.png', dpi=300, bbox_inches='tight')
print(f"  Saved: {OUT_DIR}/chromhmm_enrichment_heatmap.png")

# ===========================
# Figure 2: Bar plot comparison (GoPeaks vs MACS2)
# ===========================
print("\nGenerating bar plot comparison...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for idx, mark in enumerate(['H3K27ac', 'H3K27me3']):
    subset = enrichment_summary[enrichment_summary['Histone_Mark'] == mark]
    
    # Prepare data
    categories = sorted(subset['State_Category'].unique())
    x = np.arange(len(categories))
    width = 0.35
    
    gopeaks_scores = []
    macs2_scores = []
    
    for cat in categories:
        gopeaks_val = subset[(subset['Tool'] == 'GoPeaks') & (subset['State_Category'] == cat)]['Enrichment_Score'].values
        macs2_val = subset[(subset['Tool'] == 'MACS2') & (subset['State_Category'] == cat)]['Enrichment_Score'].values
        
        gopeaks_scores.append(gopeaks_val[0] if len(gopeaks_val) > 0 else 0)
        macs2_scores.append(macs2_val[0] if len(macs2_val) > 0 else 0)
    
    # Create bars
    axes[idx].bar(x - width/2, gopeaks_scores, width, label='GoPeaks', alpha=0.8)
    axes[idx].bar(x + width/2, macs2_scores, width, label='MACS2', alpha=0.8)
    
    axes[idx].axhline(y=1.0, color='red', linestyle='--', linewidth=1, alpha=0.5, label='Expected (no enrichment)')
    
    axes[idx].set_title(f'{mark} - Tool Comparison', fontsize=14, fontweight='bold')
    axes[idx].set_xlabel('Chromatin State Category', fontsize=12)
    axes[idx].set_ylabel('Enrichment Score', fontsize=12)
    axes[idx].set_xticks(x)
    axes[idx].set_xticklabels(categories, rotation=45, ha='right')
    axes[idx].legend()
    axes[idx].grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(f'{OUT_DIR}/chromhmm_tool_comparison.png', dpi=300, bbox_inches='tight')
print(f"  Saved: {OUT_DIR}/chromhmm_tool_comparison.png")

# ===========================
# Save summary statistics
# ===========================
print("\nSaving summary statistics...")

summary_file = f'{OUT_DIR}/chromhmm_enrichment_summary.csv'
enrichment_summary.to_csv(summary_file, index=False)
print(f"  Saved: {summary_file}")

# Key findings
print("\n" + "=" * 80)
print("KEY FINDINGS")
print("=" * 80)

for mark in ['H3K27ac', 'H3K27me3']:
    print(f"\n{mark}:")
    subset = enrichment_summary[enrichment_summary['Histone_Mark'] == mark]
    
    # Find highest enriched states for each tool
    for tool in ['GoPeaks', 'MACS2']:
        tool_data = subset[subset['Tool'] == tool].sort_values('Enrichment_Score', ascending=False)
        if len(tool_data) > 0:
            top_state = tool_data.iloc[0]
            print(f"  {tool}: Highest enrichment in '{top_state['State_Category']}' ({top_state['Enrichment_Score']:.2f}x)")

print("ChromHMM visualization complete!")

