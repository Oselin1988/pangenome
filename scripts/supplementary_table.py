#!/usr/bin/env python3
"""
Generate Supplementary Table S1: Top 50 genes with cluster, prevalence,
representative genome, and isolation source.
"""

import pandas as pd

origins = pd.read_csv('/home/lindaosaghale/eggnog_analysis/annotations/gene_origins_with_source.csv')
sel = pd.read_csv('/home/lindaosaghale/eggnog_analysis/annotations/selected_genes_xgboost.csv')

table = origins.merge(sel, on='gene', how='left')
table = table[['gene', 'importance', 'cluster', 'num_genomes',
               'representative_genome', 'isolation_source']]
table = table.sort_values(['cluster', 'importance'], ascending=[True, False])
table.columns = ['Gene', 'Importance', 'Cluster', 'Genomes carrying gene',
                 'Representative genome', 'Isolation source']

table.to_csv('/home/lindaosaghale/eggnog_analysis/annotations/Supplementary_Table_S1_50_genes.csv', index=False)
print(f"✅ Saved {len(table)} genes")
print(table.head(10).to_string(index=False))
