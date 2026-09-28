#!/usr/bin/env python3
"""
Generate additional figures: frequency vs importance and source per cluster.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

origins = pd.read_csv('/home/lindaosaghale/eggnog_analysis/annotations/gene_origins_with_source.csv')
sel = pd.read_csv('/home/lindaosaghale/eggnog_analysis/annotations/selected_genes_xgboost.csv')
df = origins.merge(sel, on='gene', how='left')

df['frequency'] = df['num_genomes'] / 1254 * 100

fig, ax = plt.subplots(figsize=(12, 8))
top20 = df.nlargest(20, 'importance')['gene'].tolist()
df['highlight'] = df['gene'].apply(lambda x: 'top20' if x in top20 else 'other')
for group, color in {'top20': 'red', 'other': 'steelblue'}.items():
    subset = df[df['highlight'] == group]
    ax.scatter(subset['frequency'], subset['importance'],
               c=color, label=group.upper(), s=50, alpha=0.7, edgecolor='black')
ax.set_xlabel('Gene frequency (%)', fontsize=14)
ax.set_ylabel('XGBoost importance', fontsize=14)
ax.legend()
ax.grid(True, linestyle='--', alpha=0.4)
plt.tight_layout()
plt.savefig('frequency_vs_importance.png', dpi=300, bbox_inches='tight')

df_known = df[df['isolation_source'].isin(['clinical', 'non_clinical'])].copy()
pivot = pd.crosstab(df_known['cluster'], df_known['isolation_source'])
all_clusters = sorted(df['cluster'].unique())
pivot = pivot.reindex(all_clusters, fill_value=0)

fig, ax = plt.subplots(figsize=(10, 6))
pivot.plot(kind='bar', stacked=True, color=['#e74c3c', '#2ecc71'],
           ax=ax, edgecolor='black')
ax.set_xlabel('Cluster', fontsize=14)
ax.set_ylabel('Number of genes', fontsize=14)
ax.tick_params(rotation=0)
plt.tight_layout()
plt.savefig('source_per_cluster.png', dpi=300, bbox_inches='tight')
print("✅ Figures saved.")
