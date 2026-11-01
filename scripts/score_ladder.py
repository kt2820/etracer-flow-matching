"""Score the conditioning ladder under the reformulated targets.

Targets (BL -> EE, one row per clone-family unit with >= 20 BL cells):
  survive         clone detected in either EE passage of its family (primary)
  growth|survived log2 growth (passages averaged), survivors only (secondary)
  growth (mixed)  log2 growth with extinction floored at 0.5 cells (the Oct 31 target)

Features, all from baseline cells only, none touching fate:
  state5   the five annotated state fractions
  mean     pseudobulk mean of per-cell log-CPM over a 400-gene panel
  kme      random-Fourier kernel mean embedding of the clone's cell distribution

The panel is fate-blind: genes with >= 2000 pooled UMI, ranked by variance of
log-CPM across the 15 tumour pseudobulks. Split: grouped 5-fold CV by clone,
20 repeats, ridge penalty by inner grouped CV, permutation null per cell.

    python3 scripts/score_ladder.py      # ~5 min; writes derived/ladder_*.csv
"""
import os, sys, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import numpy as np, pandas as pd, h5py
import scipy.sparse as sp
from scipy.stats import wilcoxon

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
import etracer as et
from reconcile_clone_tables import targets

P = et.Paths(ROOT)
PANEL_SIZE, RFF_DIM, SEED, MIN_BL = 400, 512, 7, 20
ALPHAS = (1, 10, 100, 1e3, 1e4, 1e5)
t0 = time.time()
log = lambda *a: print(f"[{time.time()-t0:6.0f}s]", *a, flush=True)

# ------------------------------------------------------------ units + targets
meta = et.load_metadata(P)
man = et.clone_manifest(meta)
D = targets(man, min_bl=MIN_BL)
D = D.rename(columns={"avg_half": "growth_mixed"})
log(f"units {len(D)}, clones {D.clone.nunique()}, survived {D.survive.mean():.1%}")

key = meta["clone"].astype(str) + "|" + meta["family"].astype(str)
unit_key = (D.clone + "|" + D.family).to_numpy()
need = (meta["stage"] == "BL").to_numpy() & meta["called"].to_numpy() & key.isin(unit_key).to_numpy()
log(f"baseline cells in units: {need.sum()}")

# ------------------------------------------- one streaming pass over the matrix
cache = P.derived / "ladder_cells.npz"
if cache.exists():
    z = np.load(cache, allow_pickle=True)
    pb, need_idx = z["pb"], z["need_idx"]
    Xneed = sp.csr_matrix((z["d"], z["i"], z["p"]), shape=tuple(z["shape"]))
    log("loaded cached pass")
else:
    with h5py.File(P.counts_h5ad, "r") as f:
        assert (f["obs/_index"].asstr()[:] == meta["cell_id"].to_numpy()).all(), "cell order differs"
        indptr = f["X/indptr"][:].astype(np.int64)
        ncell, ngene = f["X"].attrs["shape"]
        samples = meta["sample"].to_numpy(); su = np.unique(samples)
        sidx = np.searchsorted(su, samples)
        pb = np.zeros((len(su), ngene), np.float64)
        keep = []
        step = 3000
        for a in range(0, ncell, step):
            b = min(a + step, ncell)
            lo, hi = indptr[a], indptr[b]
            M = sp.csr_matrix((f["X/data"][lo:hi], f["X/indices"][lo:hi], indptr[a:b+1] - lo),
                              shape=(b - a, ngene))
            S = sp.csr_matrix((np.ones(b - a), (sidx[a:b], np.arange(b - a))), shape=(len(su), b - a))
            pb += (S @ M).toarray()
            if need[a:b].any():
                keep.append(M[need[a:b]])
            if a % 30000 == 0: log(f"  rows {a}/{ncell}")
    Xneed = sp.vstack(keep).tocsr(); need_idx = np.where(need)[0]
    np.savez_compressed(cache, pb=pb, need_idx=need_idx, d=Xneed.data, i=Xneed.indices,
                        p=Xneed.indptr, shape=np.array(Xneed.shape))
    log(f"pass done: pseudobulk {pb.shape}, baseline cells {Xneed.shape}")

# ---------------------------------------------------------- fate-blind panel
pooled = pb.sum(0)
lg = np.log2(pb / pb.sum(1, keepdims=True) * 1e6 + 1)
cand = np.where(pooled >= 2000)[0]
panel = np.sort(cand[np.argsort(-lg[:, cand].var(0))[:PANEL_SIZE]])
lib = meta["n_umi"].to_numpy()[need_idx].astype(np.float64)
Lc = np.log2(Xneed[:, panel].toarray() / lib[:, None] * 1e6 + 1).astype(np.float32)
log(f"panel {len(panel)} genes; per-cell matrix {Lc.shape}")

# ---------------------------------------------------------- per-unit features
cell_key = key.to_numpy()[need_idx]
rows = {k: np.where(cell_key == k)[0] for k in unit_key}
rng = np.random.default_rng(SEED)
mu, sg = Lc.mean(0), Lc.std(0) + 1e-6
Z = (Lc - mu) / sg
a, b = rng.choice(len(Z), 2000), rng.choice(len(Z), 2000)
gamma = 1.0 / np.median(np.linalg.norm(Z[a] - Z[b], axis=1)) ** 2
W = rng.normal(0, np.sqrt(2 * gamma), size=(Z.shape[1], RFF_DIM)); bb = rng.uniform(0, 2 * np.pi, RFF_DIM)
FEAT = {
    "state fractions (5)": D[[f"state_{s}" for s in et.STATES]].to_numpy(),
    "pseudobulk mean (400 genes)": np.stack([Lc[rows[k]].mean(0) for k in unit_key]),
    "kernel mean embedding": np.stack([np.cos(Z[rows[k]] @ W + bb).mean(0) for k in unit_key]) * np.sqrt(2 / RFF_DIM),
}
log("features:", {k: v.shape for k, v in FEAT.items()})

# ---------------------------------------------------------------- scoring
def oof(X, y, g, seed):
    rng = np.random.default_rng(seed); pred = np.full(len(y), np.nan)
    for tr, te in et.grouped_folds(g, 5, rng):
        best_a, best_s = ALPHAS[0], -np.inf
        for al in ALPHAS:
            ip = np.full(int(tr.sum()), np.nan)
            for itr, ite in et.grouped_folds(g[tr], 4, rng):
                ip[ite] = et.ridge_predict(X[tr][itr], y[tr][itr], X[tr][ite], al)
            s = et.pearson(ip, y[tr])
            if np.isfinite(s) and s > best_s: best_s, best_a = s, al
        pred[te] = et.ridge_predict(X[tr], y[tr], X[te], best_a)
    return pred

def auc(score, y):
    r = pd.Series(score).rank().to_numpy(); n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

TARGETS = {"survival (primary)": (D.survive.to_numpy(), np.ones(len(D), bool)),
           "growth | survived": (D.growth_mixed.to_numpy(), D.survive.to_numpy() == 1),
           "growth, extinction at 0.5 cells": (D.growth_mixed.to_numpy(), np.ones(len(D), bool))}
def _job(args):
    tname, fname = args
    yall, mask = TARGETS[tname]
    y, g, X = yall[mask], D.clone.to_numpy()[mask], FEAT[fname][mask]
    reps = [oof(X, y, g, s) for s in range(20)]
    r = np.array([et.pearson(p, y) for p in reps])
    nul = []
    for i in range(20):
        yp = np.random.default_rng(1000 + i).permutation(y)
        nul.append(et.pearson(oof(X, yp, g, i), yp))
    nul = np.array(nul)
    row = dict(target=tname, features=fname, n_units=int(mask.sum()), n_clones=len(np.unique(g)),
               r_mean=r.mean(), r_sd=r.std(), null_mean=nul.mean(), null_sd=nul.std())
    if tname.startswith("survival"):
        a_ = np.array([auc(p, y) for p in reps]); row.update(auc_mean=a_.mean(), auc_sd=a_.std())
    return row, r

import multiprocessing as mp
jobs = [(t, f) for t in TARGETS for f in FEAT]
jobs.sort(key=lambda j: -FEAT[j[1]].shape[1])          # biggest first
with mp.get_context("fork").Pool(4) as pool:
    res = pool.map(_job, jobs, chunksize=1)
out, per_rep = [], {}
for (t, f), (row, r) in zip(jobs, res):
    out.append(row); per_rep[(t, f)] = r
    log(t, "|", f, f"r = {row['r_mean']:.3f} +/- {row['r_sd']:.3f}  null {row['null_mean']:.3f}")
order = {k: i for i, k in enumerate(TARGETS)}; forder = {k: i for i, k in enumerate(FEAT)}
out.sort(key=lambda d: (order[d["target"]], forder[d["features"]]))
S = pd.DataFrame(out)

# ------------------------------------------------ paired: distribution vs mean
pairs = []
for tname in TARGETS:
    d = per_rep[(tname, "kernel mean embedding")] - per_rep[(tname, "pseudobulk mean (400 genes)")]
    pairs.append(dict(target=tname, kme_minus_mean=d.mean(), sd=d.std(), wins=int((d > 0).sum()),
                      of=len(d), wilcoxon_p=wilcoxon(d).pvalue))
Pq = pd.DataFrame(pairs)

# ------------------------------- reliability ceiling: same clone, other family
ceil = []
for tname, col, mask in [("survival (primary)", "survive", np.ones(len(D), bool)),
                         ("growth | survived", "growth_mixed", D.survive.to_numpy() == 1),
                         ("growth, extinction at 0.5 cells", "growth_mixed", np.ones(len(D), bool))]:
    w = D[mask].pivot_table(index="clone", columns="family", values=col).to_numpy()
    R = et.icc1([v[np.isfinite(v)] for v in w])["icc"]
    lo, hi = et.bootstrap_ci(lambda i: et.icc1([v[np.isfinite(v)] for v in w[i]])["icc"], len(w), n_boot=1000)
    ceil.append(dict(target=tname, icc=R, icc_lo=lo, icc_hi=hi, ceiling=np.sqrt(max(R, 0)),
                     ceiling_lo=np.sqrt(max(lo, 0)), ceiling_hi=np.sqrt(max(hi, 0))))
C = pd.DataFrame(ceil)

S.to_csv(P.derived / "ladder_scores.csv", index=False)
Pq.to_csv(P.derived / "ladder_paired.csv", index=False)
C.to_csv(P.derived / "ladder_ceiling.csv", index=False)
fmt = lambda x: f"{x:.3f}"
print(S.to_string(index=False, float_format=fmt)); print(Pq.to_string(index=False, float_format=fmt))
print(C.to_string(index=False, float_format=fmt)); log("done")
