#!/usr/bin/env python3
"""
Build gene origin table using Rtab only.
"""

import pandas as pd
import os

selected_file = os.path.expanduser('~/eggnog_analysis/annotations/selected_genes_xgboost.csv')
rtab_file = '/mnt/data/panaroo_merged/gene_presence_absence.Rtab'
output_file = os.path.expanduser('~/eggnog_analysis/annotations/gene_origins_table.csv')

sel_df = pd.read_csv(selected_file)
top_genes = sel_df['gene'].head(50).tolist()

rtab = pd.read_csv(rtab_file, sep='\t', index_col=0)
print(f"Rtab shape: {rtab.shape}")

results = []
for gene in top_genes:
    if gene not in rtab.index:
        results.append({'gene': gene, 'num_genomes': 0, 'representative_genome': 'Not found'})
        continue
    row = rtab.loc[gene]
    num = (row == 1).sum()
    first_genome = row[row == 1].index[0] if num > 0 else 'Not found'
    results.append({'gene': gene, 'num_genomes': num, 'representative_genome': first_genome})

pd.DataFrame(results).to_csv(output_file, index=False)
print(f"✅ Saved to {output_file}")
