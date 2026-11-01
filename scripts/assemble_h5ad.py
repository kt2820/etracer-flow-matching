#!/usr/bin/env python3
"""Turn the MatrixMarket output of convert_seurat_to_h5ad.R into an h5ad.

    python3 assemble_h5ad.py derived/matrix derived/etracer.h5ad

Needs scipy + anndata, neither of which is installed in the Cowork VM:
    pip install --user anndata scipy
Run it wherever the conversion ran (Greene), then copy only the h5ad back.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

src = Path(sys.argv[1] if len(sys.argv) > 1 else "derived/matrix")
dst = Path(sys.argv[2] if len(sys.argv) > 2 else "derived/etracer.h5ad")

from scipy.io import mmread          # noqa: E402
import anndata as ad                 # noqa: E402

genes = (src / "genes.txt").read_text().split("\n")
genes = [g for g in genes if g]
cells = (src / "cells.txt").read_text().split("\n")
cells = [c for c in cells if c]
print(f"{len(genes)} genes x {len(cells)} cells")

X = mmread(src / "counts.mtx").T.tocsr()     # cells x genes
print("matrix:", X.shape, X.dtype, f"{X.nnz} nonzero")
assert X.shape == (len(cells), len(genes))

a = ad.AnnData(X=X, obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))

# Re-verify integrality and totals before trusting it.
assert np.allclose(a.X.data, np.round(a.X.data)), "counts are not integers"
tot = pd.read_csv(src / "cell_totals.csv").set_index("cell")["total"]
mine = pd.Series(np.asarray(a.X.sum(1)).ravel(), index=a.obs_names)
both = mine.index.intersection(tot.index)
bad = int((mine.loc[both] != tot.loc[both]).sum())
print(f"per-cell totals mismatched: {bad} / {len(both)}")

a.write_h5ad(dst, compression="gzip")
print("wrote", dst)
