"""Shared helpers for the clone-conditioned flow-matching project on eTRACER.

Deliberately dependency-light: numpy + pandas only. Nothing here touches the
expression matrix, so notebooks 00-02 run on the 19 MB metadata CSV alone.

Stage names: the deposit calls them T1/T2/T3; the paper calls them
Baseline / Early-Evasion / Late-Evasion. Mapped to BL/EE/LE throughout.

Structure of the experiment, read off the sample names:
    3 BL tumours, one per founder family (a3026, a3028, a3030)
    6 EE tumours, two independent passages per family
    6 LE tumours, one per EE passage
So each family gives TWO observations of the same BL -> EE transition.
"""
from __future__ import annotations

import pickle
import re
from pathlib import Path

import numpy as np
import pandas as pd

STAGE = {"T1": "BL", "T2": "EE", "T3": "LE"}
STATES = ["1", "2", "3", "4", "5"]
DATA_DIR = "Tumor_Immune_Evasion_eTRACER"


# --------------------------------------------------------------------- paths
def find_root(start=None) -> Path:
    p = Path(start or Path.cwd()).resolve()
    for cand in [p, *p.parents]:
        if (cand / DATA_DIR).is_dir():
            return cand
    raise FileNotFoundError(f"no {DATA_DIR}/ at or above {p}")


class Paths:
    def __init__(self, root=None):
        self.root = find_root(root)
        self.data = self.root / DATA_DIR
        self.meta_csv = self.data / "scRNA" / "scRNA.metadata.csv"
        self.cell_lineage_csv = self.data / "LineageTracing" / "single_cell" / "sc_cell_lineage.csv"
        self.clone_bc_pickle = self.data / "LineageTracing" / "single_cell" / "clone_static_barcode.pickle"
        self.umi_rds = self.data / "scRNA" / "sc_UMI_count.rds"
        self.counts_h5ad = self.root / "derived" / "etracer_counts.h5ad"
        self.derived = self.root / "derived"
        self.derived.mkdir(exist_ok=True)

    def __repr__(self):
        return f"Paths(root={self.root})"


# ---------------------------------------------------------------- loading
def load_metadata(paths: Paths) -> pd.DataFrame:
    """One row per cell. `clone` is NaN where no lineage barcode was called.

    Two cell identifiers, and they are NOT interchangeable:
      `cell`     the metadata's own `cellName` (T1_a3026_AAACCCAAGACTAGAT)
      `cell_id`  the unnamed index column, which is what the count matrix
                 uses as its obs_names (T1_T1_a3026_AAACCCAAGACTAGAT-1)
    Join anything from derived/etracer_counts.h5ad on `cell_id`.
    """
    usecols = ["Unnamed: 0", "cellName", "orig.ident", "time", "lineageGrp",
               "tumor_state", "nCount_RNA", "nFeature_RNA", "harmony_clusters"]
    df = pd.read_csv(paths.meta_csv, usecols=usecols, low_memory=False)
    df = df.rename(columns={"Unnamed: 0": "cell_id",
                            "cellName": "cell", "orig.ident": "sample",
                            "lineageGrp": "clone", "tumor_state": "state",
                            "nCount_RNA": "n_umi", "nFeature_RNA": "n_genes"})
    df["stage"] = df["time"].map(STAGE)
    df["family"] = df["sample"].str.extract(r"(a30\d\d)")
    df["branch"] = [s.split(f)[-1].strip("_") or "0"
                    for s, f in zip(df["sample"], df["family"])]
    df["clone"] = df["clone"].replace("None", np.nan)
    df["state"] = df["state"].astype(str)
    df["called"] = df["clone"].notna()
    return df.drop(columns=["time"])


def load_clone_barcodes(paths: Paths) -> dict:
    with open(paths.clone_bc_pickle, "rb") as fh:
        return pickle.load(fh)


# ---------------------------------------------------------------- manifest
def clone_manifest(meta: pd.DataFrame) -> pd.DataFrame:
    """One row per (clone, sample). clone_fraction is within called cells
    of that sample, which is the only denominator the deposit supports."""
    called = meta[meta["called"]]
    denom = called.groupby("sample").size().rename("n_called_sample")

    g = called.groupby(["clone", "sample"])
    out = g.agg(n_cells=("cell", "size"), umi_median=("n_umi", "median")).reset_index()

    key = called.groupby(["clone", "sample", "state"]).size().unstack("state").fillna(0)
    key = key.reindex(columns=STATES, fill_value=0.0)
    key = key.div(key.sum(axis=1), axis=0).add_prefix("state_")
    out = out.merge(key.reset_index(), on=["clone", "sample"], how="left")

    meta_s = meta[["sample", "stage", "family", "branch"]].drop_duplicates()
    out = out.merge(meta_s, on="sample", how="left").merge(denom, on="sample", how="left")
    out["clone_fraction"] = out["n_cells"] / out["n_called_sample"]
    cols = ["clone", "sample", "stage", "family", "branch", "n_cells",
            "n_called_sample", "clone_fraction", "umi_median"] + [f"state_{s}" for s in STATES]
    return out[cols].sort_values(["family", "stage", "clone"]).reset_index(drop=True)


def sample_summary(manifest: pd.DataFrame) -> pd.DataFrame:
    """Per-sample clone diversity: how hard the bottleneck bites."""
    rows = []
    for s, d in manifest.groupby("sample"):
        p = d["clone_fraction"].to_numpy()
        p = p / p.sum()
        rows.append(dict(
            sample=s, stage=d["stage"].iloc[0], family=d["family"].iloc[0],
            n_clones=len(p), n_called=int(d["n_called_sample"].iloc[0]),
            top_clone_share=p.max(),
            eff_diversity_shannon=float(np.exp(-(p * np.log(p)).sum())),
            eff_diversity_simpson=float(1.0 / (p ** 2).sum()),
            median_cells_per_clone=float(d["n_cells"].median()),
        ))
    return pd.DataFrame(rows).sort_values(["stage", "family"]).reset_index(drop=True)


# ------------------------------------------------------------------ growth
def growth_table(manifest: pd.DataFrame, src="BL", dst="EE",
                 extinct="pseudo", pseudo_cells=0.5) -> pd.DataFrame:
    """log2 change in clone fraction from `src` to `dst`, one row per
    (clone, family, destination sample).

    A clone present at `src` in a family but absent from a `dst` sample of
    that family is the interesting case:
        extinct="pseudo"  -> observed zero, floored at pseudo_cells/n_called
        extinct="censor"  -> NaN, i.e. 'not observed' rather than 'died'
    Those two choices are not interchangeable; notebook 02 reports both.
    """
    if extinct not in {"pseudo", "censor"}:
        raise ValueError("extinct must be 'pseudo' or 'censor'")

    src_f = {(r.family, r.clone): (r.clone_fraction, r.n_cells)
             for r in manifest[manifest.stage == src].itertuples()}
    dst_m = manifest[manifest.stage == dst]
    dst_f = {(r.sample, r.clone): r.clone_fraction for r in dst_m.itertuples()}
    denom = dst_m[["sample", "family", "branch", "n_called_sample"]].drop_duplicates()

    rows = []
    for r in denom.itertuples():
        for (fam, clone), (f0, n0) in src_f.items():
            if fam != r.family:
                continue
            f1 = dst_f.get((r.sample, clone), 0.0)
            if f1 == 0.0:
                if extinct == "censor":
                    lg = np.nan
                else:
                    lg = np.log2((pseudo_cells / r.n_called_sample) / f0)
            else:
                lg = np.log2(f1 / f0)
            rows.append(dict(clone=clone, family=fam, dst_sample=r.sample,
                             branch=r.branch, f_src=f0, n_cells_src=n0,
                             f_dst=f1, detected_dst=f1 > 0, log2_growth=lg))
    return pd.DataFrame(rows)


def state_features(manifest: pd.DataFrame, stage="BL") -> pd.DataFrame:
    """Rung 1: the five state fractions per (clone, family) at `stage`."""
    d = manifest[manifest.stage == stage]
    cols = [f"state_{s}" for s in STATES]
    return d[["clone", "family", "n_cells"] + cols].copy()


# ------------------------------------------------------- stats, by hand
# scipy / sklearn are not installed in the Cowork VM, so these are closed-form.
def pearson(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return np.nan
    a, b = a[ok] - a[ok].mean(), b[ok] - b[ok].mean()
    d = np.sqrt((a ** 2).sum() * (b ** 2).sum())
    return float((a * b).sum() / d) if d > 0 else np.nan


def spearman(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return np.nan
    return pearson(pd.Series(a[ok]).rank().to_numpy(), pd.Series(b[ok]).rank().to_numpy())


def icc1(groups) -> dict:
    """One-way random-effects ICC(1): the share of variance that is between
    clones rather than between replicate measurements of the same clone.
    `groups` is an iterable of arrays, one per clone. Handles unbalanced
    group sizes with the usual k0 correction."""
    gs = [np.asarray(g, float) for g in groups]
    gs = [g[np.isfinite(g)] for g in gs]
    gs = [g for g in gs if len(g) >= 2]
    m = len(gs)
    if m < 2:
        return dict(icc=np.nan, n_groups=m, n_obs=0)
    n_i = np.array([len(g) for g in gs], float)
    N = n_i.sum()
    allv = np.concatenate(gs)
    grand = allv.mean()
    means = np.array([g.mean() for g in gs])
    msb = (n_i * (means - grand) ** 2).sum() / (m - 1)
    msw = sum(((g - g.mean()) ** 2).sum() for g in gs) / (N - m)
    k0 = (N - (n_i ** 2).sum() / N) / (m - 1)
    icc = (msb - msw) / (msb + (k0 - 1) * msw) if (msb + (k0 - 1) * msw) > 0 else np.nan
    return dict(icc=float(icc), msb=float(msb), msw=float(msw), k0=float(k0),
                n_groups=m, n_obs=int(N))


def attenuation_correct(r_obs, reliability):
    """Spearman's correction: the correlation you would see against a
    noise-free target. Only meaningful for 0 < reliability <= 1."""
    if not np.isfinite(reliability) or reliability <= 0:
        return np.nan
    return float(r_obs / np.sqrt(reliability))


def bootstrap_ci(fn, n_items, n_boot=2000, seed=0, alpha=0.05):
    """Percentile CI by resampling item indices (clones, not observations)."""
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n_items, n_items)
        v = fn(idx)
        if np.isfinite(v):
            vals.append(v)
    if not vals:
        return (np.nan, np.nan)
    return tuple(np.quantile(vals, [alpha / 2, 1 - alpha / 2]))


# ------------------------------------------------- ridge + grouped CV
def _standardize(X):
    mu, sd = X.mean(0), X.std(0)
    sd = np.where(sd == 0, 1.0, sd)
    return (X - mu) / sd, mu, sd


def ridge_predict(Xtr, ytr, Xte, alpha):
    Xs, mu, sd = _standardize(np.asarray(Xtr, float))
    ym = float(np.mean(ytr))
    p = Xs.shape[1]
    w = np.linalg.solve(Xs.T @ Xs + alpha * np.eye(p), Xs.T @ (np.asarray(ytr, float) - ym))
    return ((np.asarray(Xte, float) - mu) / sd) @ w + ym


def grouped_folds(groups, n_splits, rng):
    """Yield boolean (train, test) masks with no group on both sides."""
    groups = np.asarray(groups)
    uq = np.unique(groups)
    for chunk in np.array_split(rng.permutation(uq), n_splits):
        te = np.isin(groups, chunk)
        if te.all() or not te.any():
            continue
        yield ~te, te


def grouped_cv_r(X, y, groups, alphas=(0.1, 1, 10, 100, 1000), n_splits=5,
                 n_repeats=20, seed=0, inner_splits=4):
    """Repeated grouped CV. Returns one pooled out-of-fold Pearson r per
    repeat, so you can report mean +/- sd. Alpha chosen by inner grouped CV
    on the training folds only. Single splits on this data are unstable --
    always report the spread."""
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    groups = np.asarray(groups)
    ok = np.isfinite(y) & np.isfinite(X).all(1)
    X, y, groups = X[ok], y[ok], groups[ok]
    scores = []
    for rep in range(n_repeats):
        rng = np.random.default_rng(seed + rep)
        pred = np.full(len(y), np.nan)
        for tr, te in grouped_folds(groups, n_splits, rng):
            best_a, best_s = alphas[0], -np.inf
            for a in alphas:
                ip = np.full(int(tr.sum()), np.nan)
                for itr, ite in grouped_folds(groups[tr], inner_splits, rng):
                    ip[ite] = ridge_predict(X[tr][itr], y[tr][itr], X[tr][ite], a)
                s = pearson(ip, y[tr])
                if np.isfinite(s) and s > best_s:
                    best_s, best_a = s, a
            pred[te] = ridge_predict(X[tr], y[tr], X[te], best_a)
        scores.append(pearson(pred, y))
    return np.array(scores, float)


def permutation_null_r(X, y, groups, n_perm=20, **kw):
    """Shuffle y across groups, keeping the CV machinery identical."""
    out = []
    for i in range(n_perm):
        rng = np.random.default_rng(10_000 + i)
        out.append(grouped_cv_r(X, rng.permutation(y), groups, n_repeats=1, seed=i, **kw)[0])
    return np.array(out, float)


def report(name, scores, null=None) -> dict:
    d = dict(predictor=name, r_mean=float(np.nanmean(scores)),
             r_sd=float(np.nanstd(scores)), n_repeats=len(scores))
    if null is not None:
        d["null_mean"] = float(np.nanmean(null))
        d["null_sd"] = float(np.nanstd(null))
    return d
