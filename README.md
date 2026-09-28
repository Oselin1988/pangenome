# Pangenome-Guided Machine Learning for High-Resolution Typing and Ecological Inference of *Staphylococcus aureus*

This repository contains the data, scripts, and figures for the study:

> **Osaghale L, Ajoseh SO, Beshiru A.** *Pangenome-Guided Machine Learning for High-Resolution Typing and Ecological Inference of Staphylococcus aureus*. (Manuscript submitted to **Archives of Microbiology**).

---

## Study Overview

We analysed 1,254 high-quality *S. aureus* genomes and constructed a pangenome of 5,991 gene families using Panaroo. Functional annotation with eggNOG-mapper annotated 85.3% of clusters. Hierarchical clustering identified six natural clusters. An XGBoost classifier identified a minimal 20-gene panel achieving 0.95 test accuracy and 0.94 cross-validation accuracy (micro-average AUC = 0.99). Cluster 0 was significantly enriched for clinical isolates (15/15, p=0.0001), while Cluster 3 was exclusively non-clinical (7/7, p=0.000003).

### Key Findings

| Metric | Value |
|--------|-------|
| Genomes analysed | 1,254 |
| Gene families | 5,991 |
| Pangenome | Open (36.0% core, 17.3% accessory, 46.7% rare) |
| Clusters | 6 |
| Minimal panel | 20 genes |
| Test accuracy | 0.95 |
| CV accuracy | 0.94 (±0.03) |
| Micro-average AUC | 0.99 |
| Clinical cluster | Cluster 0 (p = 0.0001) |
| Non-clinical cluster | Cluster 3 (p = 0.000003) |

**Top discriminatory genes:**
- `group_3176` — recombinase/resolvase (Cluster 4, clinical)
- `group_159` — lactococcin (Cluster 3, non-clinical)
- `group_1316` — HxlR regulator (Cluster 3, non-clinical)

---

## Repository Structure

```
pangenome/
├── README.md
├── data/                  # Presence/absence matrix, cluster labels, annotations
├── scripts/               # Python scripts for the full analysis
├── figures/               # Publication figures (PNG, 300 DPI)
└── tables/                # Supplementary Tables S1 and S2
```

---

## Requirements

```bash
pip install pandas numpy scikit-learn xgboost matplotlib seaborn scipy biopython
```

---

## Reproducing the Analysis

```bash
python scripts/xgboost_pipeline.py
```

---

## Contact

**Linda Osaghale**  
Department of Microbiology, University of Ibadan, Ibadan, Nigeria  
Email: lindaosaghale@gmail.com

---

## Citation

> Osaghale L, Ajoseh SO, Beshiru A. Pangenome-Guided Machine Learning for High-Resolution Typing and Ecological Inference of *Staphylococcus aureus*. *Archives of Microbiology*. (Manuscript submitted).