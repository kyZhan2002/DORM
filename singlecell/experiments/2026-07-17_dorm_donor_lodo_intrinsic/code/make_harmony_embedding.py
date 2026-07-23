#!/usr/bin/env python
"""Phase 1: RNA-only, batch-robust embedding for the intrinsic-source density ratio.

Harmony-correct the dataset's 50-dim GEX PCA with batch = donor, so cell types align across
donors while composition differences (the shift DORM needs) survive. Export obs_id + Harmony
dims (H1..H50) AND raw PCA dims (PC1..PC50, baseline) to a CSV the R code joins by obs_id.
No ADT is touched (no leakage for the ADT prediction task).
"""
import os, sys
import numpy as np
import pandas as pd
import anndata as ad
import harmonypy

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
PROC = os.path.join(REPO, "singlecell", "processed", "cite_nips_RNA_proc.h5ad")
OUT = os.path.join(REPO, "singlecell", "processed", "csv_for_R", "cite_nips_EMBED_for_R.csv")

print("Loading", PROC, "(backed, obsm only)")
a = ad.read_h5ad(PROC, backed="r")
pca = np.asarray(a.obsm["GEX_X_pca"], dtype=np.float64)          # (N, 50) raw RNA PCA
N, d = pca.shape
obs_id = a.obs_names.astype(str).to_numpy()
batch = a.obs["batch"].astype(str).to_numpy()
print(f"  N={N} cells, {d} PCs, {len(np.unique(batch))} donors/batches")

np.random.seed(0)
meta = pd.DataFrame({"batch": batch})
print("Running Harmony (batch=donor) ...")
ho = harmonypy.run_harmony(pca, meta, ["batch"])
H = np.asarray(ho.Z_corr)                                         # orientation varies by version
if H.shape == (d, N):
    H = H.T
assert H.shape == (N, d), f"Harmony output shape {H.shape} != {(N, d)}"

df = pd.DataFrame({"obs_id": obs_id})
for j in range(d):
    df[f"H{j+1}"] = H[:, j]
for j in range(d):
    df[f"PC{j+1}"] = pca[:, j]
df.to_csv(OUT, index=False, float_format="%.7g")
print("Wrote", OUT, "shape", df.shape)
