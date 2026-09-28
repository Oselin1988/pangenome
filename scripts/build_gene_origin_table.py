#!/usr/bin/env python3
"""
Find source genomes and presence counts for the top 50 discriminatory genes.
"""

import pandas as pd
import os

selected_file = os.path.expanduser('~/eggnog_analysis/annotations/selected_genes_xgboost.csv')
gene_data_file = '/mnt/data/panaroo_merged/gene_data.csv'
output_file = os.path.expanduser('~/eggnog_analysis/annotations/gene_origins_table.csv')

sel_df = pd.read_csv(selected_file)
top_genes = sel_df['gene'].head(50).tolist()
print(f"Loaded {len(top_genes)} genes.")

gene_data = pd.read_csv(gene_data_file, usecols=['gene_name', 'clustering_id', 'gff_file'])

results = []
for gene in top_genes:
    rows = gene_data[gene_data['gene_name'] == gene]
    if rows.empty:
        results.append({'gene': gene, 'num_genomes': 0,
                        'representative_genome': 'Not found', 'clustering_id': 'Not found'})
        continue
    genomes = rows['gff_file'].str.replace('.gff$', '', regex=True).str.replace('^.*/', '', regex=True).tolist()
    unique_genomes = set(genomes)
    results.append({
        'gene': gene,
        'num_genomes': len(unique_genomes),
        'representative_genome': genomes[0],
        'clustering_id': rows.iloc[0]['clustering_id']
    })

pd.DataFrame(results).to_csv(output_file, index=False)
print(f"✅ Saved to {output_file}")
