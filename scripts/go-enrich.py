#!/usr/bin/env python3
"""
GO Enrichment Analysis using Enrichr API
Runs enrichment for GoPeaks and MACS2 gene lists
"""

import requests
import json
import pandas as pd
import time
import os

# Paths
GENE_LISTS_DIR = '/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/go-enrichment/gene-lists'
OUT_DIR = '/N/project/Krolab/isabella/cutandtag-peak-caller-comparison/results/go-enrichment/results'

os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 80)
print("GO Enrichment Analysis via Enrichr")
print("=" * 80)

# Enrichr API endpoints
ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
ENRICHR_RESULT_URL = 'https://maayanlab.cloud/Enrichr/enrich'

# Gene set libraries to query
LIBRARIES = [
    'GO_Biological_Process_2023',
    'GO_Molecular_Function_2023',
    'GO_Cellular_Component_2023'
]

def submit_gene_list(genes, description):
    """Submit gene list to Enrichr and get user list ID"""
    genes_str = '\n'.join(genes)
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }
    
    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception(f'Error submitting gene list: {response.text}')
    
    data = json.loads(response.text)
    return data['userListId']

def get_enrichment_results(user_list_id, library):
    """Get enrichment results for a specific library"""
    query_string = f'?userListId={user_list_id}&backgroundType={library}'
    response = requests.get(ENRICHR_RESULT_URL + query_string)
    
    if not response.ok:
        raise Exception(f'Error getting results: {response.text}')
    
    data = json.loads(response.text)
    return data[library]

def run_enrichment(gene_file, tool, mark):
    """Run full enrichment analysis for one gene list"""
    
    print(f"\n{'=' * 60}")
    print(f"Running enrichment: {tool} - {mark}")
    print(f"{'=' * 60}")
    
    # Read gene list
    with open(gene_file, 'r') as f:
        genes = [line.strip() for line in f if line.strip()]
    
    print(f"Genes: {len(genes)}")
    
    if len(genes) < 10:
        print("WARNING: Too few genes for enrichment (< 10). Skipping.")
        return None
    
    # Submit to Enrichr
    description = f'{tool}_{mark}_peaks'
    print(f"Submitting to Enrichr...")
    user_list_id = submit_gene_list(genes, description)
    print(f"User list ID: {user_list_id}")
    
    # Get results for each library
    all_results = []
    
    for library in LIBRARIES:
        print(f"  Querying {library}...")
        time.sleep(1)  # Be nice to the API
        
        results = get_enrichment_results(user_list_id, library)
        
        # Parse results
        for result in results[:20]:  # Top 20 terms
            rank, term, p_value, z_score, combined_score, genes_overlap, adj_p_value, old_p, old_adj = result
            
            all_results.append({
                'Tool': tool,
                'Histone_Mark': mark,
                'Library': library,
                'Rank': rank,
                'Term': term,
                'P_value': p_value,
                'Adjusted_P_value': adj_p_value,
                'Combined_Score': combined_score,
                'Genes': genes_overlap
            })
    
    return pd.DataFrame(all_results)

# ===========================
# Run enrichment for all conditions
# ===========================

all_enrichment_results = []

gene_lists = {
    ('GoPeaks', 'H3K27ac'): f'{GENE_LISTS_DIR}/GoPeaks_H3K27ac_genes.txt',
    ('GoPeaks', 'H3K27me3'): f'{GENE_LISTS_DIR}/GoPeaks_H3K27me3_genes.txt',
    ('MACS2', 'H3K27ac'): f'{GENE_LISTS_DIR}/MACS2_H3K27ac_genes.txt',
    ('MACS2', 'H3K27me3'): f'{GENE_LISTS_DIR}/MACS2_H3K27me3_genes.txt',
}

for (tool, mark), gene_file in gene_lists.items():
    if not os.path.exists(gene_file):
        print(f"\nWARNING: Gene list not found: {gene_file}")
        continue
    
    try:
        results_df = run_enrichment(gene_file, tool, mark)
        if results_df is not None:
            all_enrichment_results.append(results_df)
    except Exception as e:
        print(f"ERROR running enrichment for {tool} {mark}: {e}")
        continue

# ===========================
# Combine and save results
# ===========================

if all_enrichment_results:
    combined_df = pd.concat(all_enrichment_results, ignore_index=True)
    
    output_file = f'{OUT_DIR}/go_enrichment_results.csv'
    combined_df.to_csv(output_file, index=False)
    print(f"\n✓ Results saved: {output_file}")
    
    # Summary
    print("\n" + "=" * 80)
    print("ENRICHMENT SUMMARY")
    print("=" * 80)
    
    for tool in ['GoPeaks', 'MACS2']:
        for mark in ['H3K27ac', 'H3K27me3']:
            subset = combined_df[(combined_df['Tool'] == tool) & 
                                (combined_df['Histone_Mark'] == mark) &
                                (combined_df['Library'] == 'GO_Biological_Process_2023')]
            
            if len(subset) > 0:
                print(f"\n{tool} - {mark} (Top 5 BP terms):")
                top5 = subset.nsmallest(5, 'Adjusted_P_value')
                for _, row in top5.iterrows():
                    print(f"  • {row['Term'][:60]} (p={row['Adjusted_P_value']:.2e})")
    
    print("\n" + "=" * 80)
    print("✓ GO enrichment complete!")
    print("=" * 80)
    
else:
    print("\nERROR: No enrichment results generated")
