#!/usr/bin/env python3
"""
Visualize ENCODE comparison results
Creates bar plots and summary statistics for peak caller performance
"""
 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
 
# Paths
H3K27AC_FILE = '/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-comparison/encode-overlap/encode_comparison_H3K27ac.csv'
H3K27ME3_FILE = '/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-comparison/encode-overlap/encode_comparison_H3K27me3.csv'
OUT_DIR = '/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/peak-comparison/figures'
 
import os
os.makedirs(OUT_DIR, exist_ok=True)
# Set plotting style
sns.set_style("whitegrid")
sns.set_palette("Set2")
 
print("=" * 80)
print("ENCODE Comparison Visualization")
print("=" * 80)
 
# Read data
print("\nReading ENCODE comparison results...")
h3k27ac = pd.read_csv(H3K27AC_FILE)
h3k27me3 = pd.read_csv(H3K27ME3_FILE)
 
# Combine
all_data = pd.concat([h3k27ac, h3k27me3], ignore_index=True)
 
print(f"\nTotal comparisons: {len(all_data)}")
print(f"  H3K27ac: {len(h3k27ac)}")
print(f"  H3K27me3: {len(h3k27me3)}")

# Print summary statistics
print("SUMMARY STATISTICS")
print("=" * 80)
 
for mark in ['H3K27ac', 'H3K27me3']:
    print(f"\n{mark}:")
    print("-" * 40)
    subset = all_data[all_data['Histone_Mark'] == mark]
    summary = subset.groupby('Tool')[['Recall_%', 'Precision_%']].agg(['mean', 'std'])
    print(summary)

# Create figure for % recall by tool
print("Creating figure for recall by tool...")
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
 
for idx, mark in enumerate(['H3K27ac', 'H3K27me3']):
    subset = all_data[all_data['Histone_Mark'] == mark]
    
    # Box plot
    sns.boxplot(data=subset, x='Tool', y='Recall_%', ax=axes[idx])
    axes[idx].set_title(f'{mark} - ENCODE Recall %', fontsize=14, fontweight='bold')
    axes[idx].set_xlabel('Peak Caller', fontsize=12)
    axes[idx].set_ylabel('Recall (% ENCODE peaks recovered)', fontsize=12)
    axes[idx].tick_params(axis='x', rotation=45)
    axes[idx].set_ylim(0, 100)
    
    # Add sample size
    n_samples = len(subset['Sample'].unique())
    axes[idx].text(0.02, 0.98, f'n={n_samples} samples', 
                   transform=axes[idx].transAxes, 
                   verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
plt.tight_layout()
plt.savefig(f'{OUT_DIR}/encode_recall_by_tool.png', dpi=300, bbox_inches='tight')
print(f"  Saved: {OUT_DIR}/encode_recall_by_tool.png")
 
# Create figure for % precision by tool
print("Creating figure for precision by tool...")
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
 
for idx, mark in enumerate(['H3K27ac', 'H3K27me3']):
    subset = all_data[all_data['Histone_Mark'] == mark]
    
    sns.boxplot(data=subset, x='Tool', y='Precision_%', ax=axes[idx])
    axes[idx].set_title(f'{mark} - Precision %', fontsize=14, fontweight='bold')
    axes[idx].set_xlabel('Peak Caller', fontsize=12)
    axes[idx].set_ylabel('Precision (% called peaks in ENCODE)', fontsize=12)
    axes[idx].tick_params(axis='x', rotation=45)
    axes[idx].set_ylim(0, 100)
    
    n_samples = len(subset['Sample'].unique())
    axes[idx].text(0.02, 0.98, f'n={n_samples} samples', 
                   transform=axes[idx].transAxes, 
                   verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
 
plt.tight_layout()
plt.savefig(f'{OUT_DIR}/encode_precision_by_tool.png', dpi=300, bbox_inches='tight')
print(f"  Saved: {OUT_DIR}/encode_precision_by_tool.png")

# Create scatter plot for recall vs precision
print("Creating figure for scatter plot of recal vs precision...")
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
 
tools = all_data['Tool'].unique()
colors = sns.color_palette("Set2", n_colors=len(tools))
tool_colors = dict(zip(tools, colors))
 
for idx, mark in enumerate(['H3K27ac', 'H3K27me3']):
    subset = all_data[all_data['Histone_Mark'] == mark]
    
    for tool in tools:
        tool_data = subset[subset['Tool'] == tool]
        axes[idx].scatter(tool_data['Recall_%'], tool_data['Precision_%'], 
                         label=tool, alpha=0.6, s=100, color=tool_colors[tool])
    
    axes[idx].set_title(f'{mark} - Recall vs Precision', fontsize=14, fontweight='bold')
    axes[idx].set_xlabel('Recall (% ENCODE peaks recovered)', fontsize=12)
    axes[idx].set_ylabel('Precision (% called peaks in ENCODE)', fontsize=12)
    axes[idx].legend()
    axes[idx].set_xlim(0, 100)
    axes[idx].set_ylim(0, 100)
    axes[idx].grid(True, alpha=0.3)
 
plt.tight_layout()
plt.savefig(f'{OUT_DIR}/encode_recall_vs_precision.png', dpi=300, bbox_inches='tight')
print(f"  Saved: {OUT_DIR}/encode_recall_vs_precision.png")
 
# Create figure for average performance comparison
print("Creating figure for average performance comparision...")
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
 
for idx, mark in enumerate(['H3K27ac', 'H3K27me3']):
    subset = all_data[all_data['Histone_Mark'] == mark]
    summary = subset.groupby('Tool')[['Recall_%', 'Precision_%']].mean().reset_index()
    
    x = np.arange(len(summary))
    width = 0.35
    
    axes[idx].bar(x - width/2, summary['Recall_%'], width, label='Recall %', alpha=0.8)
    axes[idx].bar(x + width/2, summary['Precision_%'], width, label='Precision %', alpha=0.8)
    
    axes[idx].set_title(f'{mark} - Average Performance', fontsize=14, fontweight='bold')
    axes[idx].set_ylabel('Percentage', fontsize=12)
    axes[idx].set_xticks(x)
    axes[idx].set_xticklabels(summary['Tool'], rotation=45, ha='right')
    axes[idx].legend()
    axes[idx].set_ylim(0, 100)
    axes[idx].grid(axis='y', alpha=0.3)
 
plt.tight_layout()
plt.savefig(f'{OUT_DIR}/encode_average_performance.png', dpi=300, bbox_inches='tight')
print(f"  Saved: {OUT_DIR}/encode_average_performance.png")

# Save summary table
summary_table = all_data.groupby(['Histone_Mark', 'Tool']).agg({
    'Recall_%': ['mean', 'std', 'min', 'max'],
    'Precision_%': ['mean', 'std', 'min', 'max'],
    'Total_Peaks': 'mean',
    'ENCODE_Overlap': 'mean'
}).round(2)
 
summary_file = f'{OUT_DIR}/encode_summary_statistics.csv'
summary_table.to_csv(summary_file)
print(f"  Saved: {summary_file}")
 

