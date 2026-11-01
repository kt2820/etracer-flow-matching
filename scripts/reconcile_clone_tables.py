"""Reconcile rung 1 (5 BL state fractions -> BL->EE clone growth, held-out
clones, grouped 5-fold CV x 20): Oct 30 r = 0.323 vs Oct 31 r = 0.263.

Two things differed between those runs:
  1. the clone table: Zenodo 16406856 (sc.metadata.BLEELE.csv) vs
     Zenodo 19230230 (scRNA.metadata.csv, the deposit in this folder);
  2. the growth definition: EE passages pooled vs averaged, and extinction
     floored at a fraction of 1e-5 vs at 0.5 cells.

Usage:
    python3 scripts/reconcile_clone_tables.py [path/to/sc.metadata.BLEELE.csv]
Without the old file, only part 2 runs. The old file can be pulled out of the
16406856 zip without downloading it all:
    pip install remotezip
    python3 -c "from remotezip import RemoteZip as Z; Z('https://zenodo.org/api/records/16406856/files/eTRACER.zip/content').extract('eTRACER/scRNA/sc.metadata.BLEELE.csv','derived/')"
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
import etracer as et

P = et.Paths(ROOT)
Xc = [f"state_{s}" for s in et.STATES]


def table_from(clones, meta):
    m = meta.copy(); m["clone"] = clones; m["called"] = m["clone"].notna(); return m


def targets(man, min_bl=20):
    """Every growth definition side by side, one row per (clone, family)."""
    bl, ee = man[man.stage == "BL"], man[man.stage == "EE"]
    rows = []
    for fam in sorted(bl.family.unique()):
        b = bl[(bl.family == fam) & (bl.n_cells >= min_bl)]
        e = ee[ee.family == fam]
        tot = e.drop_duplicates("sample").set_index("sample").n_called_sample
        for r in b.itertuples():
            per = {s: e[(e["sample"] == s) & (e.clone == r.clone)].n_cells.sum() for s in tot.index}
            pooled, f0 = sum(per.values()) / tot.sum(), r.clone_fraction
            rows.append(dict(
                clone=r.clone, family=fam,
                avg_half=np.mean([np.log2(max(per[s], 0.5) / tot[s] / f0) for s in tot.index]),   # Oct 31
                pool_half=np.log2(max(pooled, 0.5 / tot.sum()) / f0),
                avg_1e5=np.mean([np.log2((per[s] / tot[s] + 1e-5) / (f0 + 1e-5)) for s in tot.index]),
                pool_1e5=np.log2((pooled + 1e-5) / (f0 + 1e-5)),                                  # Oct 30
                survive=float(pooled > 0)))
    return pd.DataFrame(rows).merge(et.state_features(man, "BL").drop(columns="n_cells"),
                                    on=["clone", "family"])


def score(D, col, mask=None):
    d = D if mask is None else D[mask]
    s = et.grouped_cv_r(d[Xc].to_numpy(), d[col].to_numpy(), d["clone"].to_numpy(), n_repeats=20, seed=0)
    return s.mean(), s.std()


def main():
    meta = et.load_metadata(P)
    tables = {"19230230 (this folder)": meta}
    if len(sys.argv) > 1:
        old = pd.read_csv(sys.argv[1], index_col=0, low_memory=False)
        assert (old.index == meta["cell_id"].values).all(), "cell order differs"
        o, n = old["lineageGrp"].values, meta["clone"].values
        both = pd.notna(o) & pd.notna(n)
        ct = pd.crosstab(pd.Series(n[both]), pd.Series(o[both]))
        print(f"clone tables: {both.sum()} cells called in both; mapping one-to-one: "
              f"{ct.idxmax(axis=1).nunique() == len(ct)}, min purity {(ct.max(axis=1) / ct.sum(axis=1)).min():.3f}")
        for st in ["BL", "EE", "LE"]:
            m = (meta["stage"] == st).values
            print(f"  {st}: call rate {pd.notna(n[m]).mean():.3f} (19230230) vs {pd.notna(o[m]).mean():.3f} (16406856)")
        tables["16406856 (Oct 30)"] = table_from(old["lineageGrp"].values, meta)

    rows = []
    for tn, m in tables.items():
        D = targets(et.clone_manifest(m))
        ext = D.survive == 0
        for col in ["avg_half", "pool_half", "avg_1e5", "pool_1e5", "survive"]:
            r, sd = score(D, col)
            rows.append(dict(clone_table=tn, target=col, n_units=len(D), r_mean=r, r_sd=sd,
                             survival_share_of_variance=np.corrcoef(D[col], D.survive)[0, 1] ** 2))
        r, sd = score(D, "avg_half", ~ext)
        rows.append(dict(clone_table=tn, target="avg_half, survivors only", n_units=int((~ext).sum()),
                         r_mean=r, r_sd=sd, survival_share_of_variance=np.nan))
    out = pd.DataFrame(rows)
    out.to_csv(P.derived / "rung1_decomposition.csv", index=False)
    print(out.to_string(index=False, float_format=lambda x: f"{x:.3f}"))


if __name__ == "__main__":
    main()
