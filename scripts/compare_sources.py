#!/usr/bin/env python3
"""
Compare the source distribution of the 50-gene panel vs the overall dataset.
"""

import pandas as pd
from scipy.stats import chi2_contingency, binomtest

overall = {'clinical': 509, 'non_clinical': 55, 'other': 333, 'unknown': 358}
total = sum(overall.values())

observed = {'clinical': 17, 'non_clinical': 7, 'other': 2, 'unknown': 24}
n_genes = sum(observed.values())

expected = {cat: n_genes * (overall[cat] / total) for cat in overall}
obs_list = [observed[cat] for cat in overall]
exp_list = [expected[cat] for cat in overall]

chi2, p, dof, _ = chi2_contingency([obs_list, exp_list], correction=False)
print("=== Comparison of 50-gene source distribution vs overall ===")
print(f"Chi-square test: chi2 = {chi2:.3f}, p = {p:.5f}")

clin_obs = observed['clinical']
nonclin_obs = observed['non_clinical']
total_known = clin_obs + nonclin_obs
overall_clin_prop = overall['clinical'] / (overall['clinical'] + overall['non_clinical'])
p_binom = binomtest(clin_obs, total_known, overall_clin_prop, alternative='two-sided').pvalue

print(f"\nBinomial test (clinical vs non-clinical):")
print(f"  Observed clinical: {clin_obs}/{total_known} ({clin_obs/total_known*100:.1f}%)")
print(f"  Expected clinical: {overall_clin_prop*100:.1f}%")
print(f"  p = {p_binom:.5f}")
