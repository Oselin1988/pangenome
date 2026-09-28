#!/usr/bin/env python3
"""
Publication‑ready XGBoost pipeline.
Handles singleton clusters by filtering and remapping labels to consecutive integers.
"""

import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, cross_val_score, KFold, GridSearchCV
from sklearn.metrics import accuracy_score, classification_report
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import pairwise_distances, silhouette_score
import xgboost as xgb

# ---------- 1. Load data ----------
print("Loading data...")
X_df = pd.read_csv('genome_by_gene_matrix.csv', index_col=0)
print(f"Matrix: {X_df.shape[0]} genomes × {X_df.shape[1]} genes")
X_np = np.ascontiguousarray(X_df.values.astype(np.float32))
feature_names = X_df.columns.tolist()

# ---------- 2. Generate clustering labels (if missing) ----------
if not os.path.exists('typing_labels.csv'):
    print("Generating clustering labels...")
    dist = pairwise_distances(X_np, metric='jaccard')
    best_score = -1
    best_labels = None
    best_threshold = None
    for thresh in np.arange(0.15, 0.5, 0.05):
        clustering = AgglomerativeClustering(n_clusters=None, distance_threshold=thresh,
                                             metric='precomputed', linkage='average')
        labels = clustering.fit_predict(dist)
        n_clusters = len(set(labels))
        if n_clusters < 2:
            continue
        score = silhouette_score(dist, labels, metric='precomputed')
        print(f"  Threshold {thresh:.2f}: {n_clusters} clusters, silhouette={score:.3f}")
        if score > best_score:
            best_score = score
            best_labels = labels
            best_threshold = thresh
    print(f"Best threshold: {best_threshold:.2f} (silhouette={best_score:.3f}), clusters: {len(set(best_labels))}")
    pd.Series(best_labels, index=X_df.index, name='cluster').to_csv('typing_labels.csv')
    print("✅ typing_labels.csv saved")
else:
    print("typing_labels.csv exists.")

# ---------- 3. Load labels and filter singleton clusters ----------
y_df = pd.read_csv('typing_labels.csv', index_col=0)['cluster']
cluster_sizes = y_df.value_counts()
print("\nCluster sizes before filtering:")
print(cluster_sizes)

singletons = cluster_sizes[cluster_sizes < 2].index.tolist()
if singletons:
    print(f"⚠️ Removing {len(singletons)} singleton clusters: {singletons}")
    keep_mask = ~y_df.isin(singletons)
    X_df = X_df.loc[keep_mask]
    y_df = y_df.loc[keep_mask]
    print(f"Remaining genomes: {X_df.shape[0]}")
    # Remap labels to consecutive integers
    unique_labels = sorted(y_df.unique())
    label_map = {old: new for new, old in enumerate(unique_labels)}
    y_df = y_df.map(label_map)
    print(f"Labels remapped: {unique_labels} -> {sorted(y_df.unique())}")
    y_df.to_csv('typing_labels_filtered.csv')
    print("✅ Filtered and remapped labels saved to typing_labels_filtered.csv")
else:
    print("No singleton clusters found.")
    # Still ensure labels are consecutive (just in case)
    unique_labels = sorted(y_df.unique())
    if list(unique_labels) != list(range(len(unique_labels))):
        label_map = {old: new for new, old in enumerate(unique_labels)}
        y_df = y_df.map(label_map)
        print(f"Labels remapped: {unique_labels} -> {sorted(y_df.unique())}")
        y_df.to_csv('typing_labels_remapped.csv')
        print("✅ Remapped labels saved to typing_labels_remapped.csv")

# ---------- 4. Prepare data ----------
X_np = np.ascontiguousarray(X_df.values.astype(np.float32))
y_np = np.ascontiguousarray(y_df.values.astype(int))
feature_names = X_df.columns.tolist()

print(f"Final data shape: {X_np.shape}, labels: {len(np.unique(y_np))} clusters")

# ---------- 5. Split data ----------
X_train, X_test, y_train, y_test = train_test_split(
    X_np, y_np, test_size=0.2, random_state=42  # no stratification to avoid rare class issues
)

# ---------- 6. XGBoost tuning ----------
print("\n--- Tuning XGBoost ---")
param_grid = {
    'max_depth': [4, 6, 8],
    'learning_rate': [0.05, 0.1, 0.2],
    'subsample': [0.8, 1.0],
    'colsample_bytree': [0.8, 1.0],
    'n_estimators': [100, 200]
}
xgb_clf = xgb.XGBClassifier(objective='multi:softprob', eval_metric='mlogloss',
                            random_state=42, use_label_encoder=False)
cv = KFold(n_splits=5, shuffle=True, random_state=42)
grid = GridSearchCV(xgb_clf, param_grid, cv=cv, scoring='accuracy', n_jobs=-1, verbose=1)
grid.fit(X_train, y_train)
best_xgb = grid.best_estimator_
print(f"Best params: {grid.best_params_}")
print(f"CV accuracy: {grid.best_score_:.3f}")

y_pred = best_xgb.predict(X_test)
test_acc = accuracy_score(y_test, y_pred)
print(f"Test accuracy: {test_acc:.3f}")
print("\nClassification report:\n", classification_report(y_test, y_pred))

# ---------- 7. Feature selection ----------
importances = best_xgb.feature_importances_
top_n = 50
indices = np.argsort(importances)[::-1][:top_n]
selected_genes = [feature_names[i] for i in indices]

print(f"\nSelected {len(selected_genes)} genes (top {top_n} by importance)")

X_train_sel = X_train[:, indices]
X_test_sel = X_test[:, indices]
xgb_sel = xgb.XGBClassifier(**grid.best_params_, objective='multi:softprob',
                            eval_metric='mlogloss', random_state=42, use_label_encoder=False)
xgb_sel.fit(X_train_sel, y_train)
y_pred_sel = xgb_sel.predict(X_test_sel)
acc_sel = accuracy_score(y_test, y_pred_sel)
cv_sel = cross_val_score(xgb_sel, X_train_sel, y_train, cv=5, scoring='accuracy')
print(f"Selected model test accuracy: {acc_sel:.3f} (CV: {cv_sel.mean():.3f} +/- {cv_sel.std():.3f})")

# Save selected genes
sel_df = pd.DataFrame({
    'gene': selected_genes,
    'importance': xgb_sel.feature_importances_
}).sort_values('importance', ascending=False)
sel_df.to_csv('selected_genes_xgboost.csv', index=False)
print("✅ selected_genes_xgboost.csv saved")

# ---------- 8. Feature importance plot ----------
plt.figure(figsize=(14, 10))
imp = xgb_sel.feature_importances_
idx = np.argsort(imp)[::-1][:20]
plt.barh(range(len(idx)), imp[idx])
plt.yticks(range(len(idx)), [selected_genes[i] for i in idx])
plt.xlabel('Importance')
plt.title('Top 20 predictive genes (XGBoost)')
plt.tight_layout()
plt.savefig('feature_importance_xgboost.png', dpi=300)
print("✅ Feature importance plot saved to feature_importance_xgboost.png")

# ---------- 9. Annotations ----------
print("\n--- Top predictive genes with functional annotations ---")
annot = pd.read_csv('s_aureus_eggnog_annotations.csv', index_col='query')
merged = sel_df.head(10).merge(annot, left_on='gene', right_index=True, how='left')
print(merged[['gene', 'importance', 'COG_category', 'KEGG_Pathway']].to_string(index=False))

# ---------- 10. Per‑cluster enrichment ----------
print("\n--- Per-cluster enrichment of top 10 genes ---")
labels = y_df  # already filtered and remapped
top_genes = sel_df['gene'].head(10).tolist()
for cluster in sorted(set(labels)):
    c_genomes = labels[labels == cluster].index
    sub = X_df.loc[c_genomes, top_genes]
    freq = (sub == 1).mean(axis=0)
    print(f"\nCluster {cluster} (n={len(c_genomes)}):")
    print(freq.sort_values(ascending=False).to_string())

print("\n✅ Pipeline completed successfully.")
