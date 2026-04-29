#!/usr/bin/env python3
"""
GO Enrichment Visualization - GoPeaks vs MACS2
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Paths
INPUT_CSV = "/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/go-enrichment/results/go_enrichment_results.csv"
OUTPUT_PNG = "/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/go-enrichment/results/go_enrichment_comparison.png"

#Load data 
df = pd.read_csv(INPUT_CSV)

# Filter to H3K27ac, GO Biological Process only
df_bp = df[
    (df['Histone_Mark'] == 'H3K27ac') &
    (df['Library'] == 'GO_Biological_Process_2023')
].copy()

# Get top 5 terms per tool by adjusted p-value
TOP_N = 5

gopeaks_df = (df_bp[df_bp['Tool'] == 'GoPeaks']
              .sort_values('Adjusted_P_value')
              .head(TOP_N)
              .reset_index(drop=True))

macs2_df = (df_bp[df_bp['Tool'] == 'MACS2']
            .sort_values('Adjusted_P_value')
            .head(TOP_N)
            .reset_index(drop=True))

# Compute -log10 adjusted p-values
gopeaks_df['neg_log_p'] = -np.log10(gopeaks_df['Adjusted_P_value'])
macs2_df['neg_log_p']   = -np.log10(macs2_df['Adjusted_P_value'])

# Clean up term names for display (strip GO accession)
def clean_term(term):
    return term.split(' (GO:')[0].replace(' ', '\n', 1) if len(term) > 25 else term

gopeaks_df['label'] = gopeaks_df['Term'].apply(clean_term)
macs2_df['label']   = macs2_df['Term'].apply(clean_term)

# Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
fig.patch.set_facecolor('white')

GOPEAKS_COLOR = '#2E86AB'
MACS2_COLOR   = '#E84855'
SIG_LINE      = '#666666'
SIG_THRESHOLD = -np.log10(0.05)

def make_panel(ax, data, color, title):
    bars = ax.barh(range(len(data)), data['neg_log_p'],
                   color=color, alpha=0.85, edgecolor='white', linewidth=0.5)
    ax.set_yticks(range(len(data)))
    ax.set_yticklabels(data['label'], fontsize=10)
    ax.set_xlabel('-log₁₀(Adjusted P-value)', fontsize=11)
    ax.set_title(title, fontsize=12, fontweight='bold', color=color)
    ax.axvline(x=SIG_THRESHOLD, color=SIG_LINE, linestyle='--',
               linewidth=1.2, alpha=0.7, label='p = 0.05')
    ax.invert_yaxis()
    ax.set_xlim(0, max(data['neg_log_p'].max() * 1.2, SIG_THRESHOLD * 1.5))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(fontsize=9, loc='lower right')
    # Value labels
    for bar, val in zip(bars, data['neg_log_p']):
        ax.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height() / 2,
                f'{val:.1f}', va='center', fontsize=9, color='#333333')

make_panel(ax1, gopeaks_df, GOPEAKS_COLOR,
           'GoPeaks H3K27ac\nTop GO Biological Processes')
make_panel(ax2, macs2_df,   MACS2_COLOR,
           'MACS2 H3K27ac\nTop GO Biological Processes')

plt.suptitle(
    'Figure 3. GO Biological Process Enrichment: GoPeaks vs MACS2 (H3K27ac)',
    fontsize=11, fontweight='bold', y=1.01
)
plt.tight_layout()
plt.savefig(OUTPUT_PNG, dpi=200, bbox_inches='tight', facecolor='white')
print(f"Saved to {OUTPUT_PNG}")
