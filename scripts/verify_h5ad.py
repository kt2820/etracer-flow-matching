#!/usr/bin/env python3
"""Verify the converted h5ad against the deposited metadata.

The check that matters: per-cell UMI totals recomputed from the matrix must
equal nCount_RNA exactly. That simultaneously confirms the values, the column
order and the cell-name mapping. nFeature_RNA checks the index arrays.

    python3 verify_h5ad.py derived/etracer_counts.h5ad <meta.csv>
"""
import sys

import h5py
import numpy as np
import pandas as pd

h5, meta_csv = sys.argv[1], sys.argv[2]
f = h5py.File(h5, "r")
print("root attrs:", {k: f.attrs[k] for k in f.attrs})
print("X attrs:   ", {k: f["X"].attrs[k] for k in f["X"].attrs})

p = f["X/indptr"][:].astype(np.int64)
cells = f["obs/_index"].asstr()[:]
genes = f["var/_index"].asstr()[:]
n_cells, n_genes = len(p) - 1, len(genes)
nnz = int(p[-1])
print(f"\n{n_cells:,} cells x {n_genes:,} genes, nnz {nnz:,}")

data, ind = f["X/data"], f["X/indices"]
sums = np.zeros(n_cells, dtype=np.float64)
nfeat = (p[1:] - p[:-1]).astype(np.int64)
gmax = -1
BLOCK = 1 << 23
c = 0
while c < n_cells:
    c2 = int(np.searchsorted(p, p[c] + BLOCK, side="right"))
    c2 = min(max(c2, c + 1), n_cells)
    s, e = int(p[c]), int(p[c2])
    if e > s:
        blk = data[s:e]
        bounds = (p[c:c2] - s).astype(np.int64)
        seg = np.add.reduceat(blk, bounds).astype(np.float64)
        seg[nfeat[c:c2] == 0] = 0.0            # reduceat mis-handles empty spans
        sums[c:c2] = seg
        gmax = max(gmax, int(ind[s:e].max()))
    c = c2

md = pd.read_csv(meta_csv, index_col=0, low_memory=False)
print(f"metadata rows: {len(md):,}")

same_order = bool(len(md) == n_cells and (md.index.to_numpy() == cells).all())
print(f"\ncell names identical and in the same order: {same_order}")
if not same_order:
    aligned = md.reindex(cells)
    print("  -> reindexed metadata onto matrix order; "
          f"missing: {int(aligned['nCount_RNA'].isna().sum())}")
    md = aligned

bad_count = int((md["nCount_RNA"].to_numpy(float) != sums).sum())
bad_feat = int((md["nFeature_RNA"].to_numpy(np.int64) != nfeat).sum())
print(f"nCount_RNA  mismatches: {bad_count} / {n_cells}")
print(f"nFeature_RNA mismatches: {bad_feat} / {n_cells}")
print(f"max gene index {gmax} (< {n_genes}: {gmax < n_genes})")
print(f"genes unique: {len(set(genes)) == n_genes} | cells unique: {len(set(cells)) == n_cells}")
print(f"total UMI: {sums.sum():,.0f}")

try:
    import anndata as ad
    a = ad.read_h5ad(h5)
    print(f"\nanndata reads it: {a}")
    print(f"  X dtype {a.X.dtype}, format {type(a.X).__name__}")
    sub = a[:200].X.sum(1).A.ravel()
    print(f"  first-200 round-trip matches: {bool(np.allclose(sub, sums[:200]))}")
except ImportError:
    print("\n(anndata not installed; h5py-level checks only)")

ok = bad_count == 0 and bad_feat == 0 and gmax < n_genes
print("\n" + ("PASS: the matrix is faithful to the deposit"
              if ok else "FAIL: see mismatches above"))
