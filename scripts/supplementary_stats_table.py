#!/usr/bin/env python3
"""
Generate Supplementary Table S2: Summary of all statistical tests.
"""

import pandas as pd

rows = [
    {'Analysis': 'Overall source distribution',
     'Category / Test': 'Chi-square (4 categories)',
     'Comparison / Group': '50 genes vs. 1,255 genomes',
     'Statistic / Count': 'chi2 = 13.60', 'p-value': '0.0035',
     'Interpretation / Enrichment': 'Significant difference'},
    {'Analysis': 'Overall source distribution',
     'Category / Test': 'Binomial (Clinical vs Non-clinical)',
     'Comparison / Group': '50 genes (known sources)',
     'Statistic / Count': 'Observed 70.8% (17/24) vs Expected 90.2%',
     'p-value': '0.0065',
     'Interpretation / Enrichment': 'Significant depletion of clinical genes'},
    {'Analysis': 'Cluster-wise enrichment (Fisher)',
     'Category / Test': 'Cluster 0 vs. others',
     'Comparison / Group': 'Clinical (15) vs Non-clinical (0)',
     'Statistic / Count': '100% clinical', 'p-value': '0.0001',
     'Interpretation / Enrichment': 'Clinical-enriched'},
    {'Analysis': 'Cluster-wise enrichment (Fisher)',
     'Category / Test': 'Cluster 3 vs. others',
     'Comparison / Group': 'Clinical (0) vs Non-clinical (7)',
     'Statistic / Count': '100% non-clinical', 'p-value': '0.000003',
     'Interpretation / Enrichment': 'Non-clinical-enriched'},
    {'Analysis': 'Cluster-wise enrichment (Fisher)',
     'Category / Test': 'Cluster 4 vs. others',
     'Comparison / Group': 'Clinical (2) vs Non-clinical (0)',
     'Statistic / Count': '100% clinical (n=2)', 'p-value': '1.00',
     'Interpretation / Enrichment': 'Not significant'},
    {'Analysis': 'Cluster-wise enrichment (Fisher)',
     'Category / Test': 'Clusters 1, 2, 5',
     'Comparison / Group': 'All sources unknown',
     'Statistic / Count': '0 known-source genes', 'p-value': 'N/A',
     'Interpretation / Enrichment': 'Excluded from analysis'},
    {'Analysis': 'COG functional enrichment',
     'Category / Test': 'Fisher exact test (all COG categories)',
     'Comparison / Group': 'Clinical vs Non-clinical genes',
     'Statistic / Count': 'No category enriched', 'p-value': '>= 0.1 (all)',
     'Interpretation / Enrichment': 'No significant enrichment'},
    {'Analysis': 'Domain enrichment (top 20)',
     'Category / Test': 'Fisher exact test',
     'Comparison / Group': 'Clinical vs Non-clinical genes',
     'Statistic / Count': 'No domain enriched', 'p-value': '>= 0.1 (all)',
     'Interpretation / Enrichment': 'No significant enrichment'},
]

df = pd.DataFrame(rows)
df.to_csv('/home/lindaosaghale/eggnog_analysis/annotations/Supplementary_Table_Stats.csv', index=False)
print(f"✅ Saved {len(df)} rows")
print(df.to_string(index=False))
