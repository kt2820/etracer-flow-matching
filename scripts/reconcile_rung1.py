"""Which setup produced which rung-1 number? Rung 1 (5 BL state fractions ->
BL->EE clone log2 growth), held-out-clone grouped 5-fold CV x 20, every
combination of the three choices that differed between sessions."""
import sys, itertools
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import etracer as et

P = et.Paths(Path(__file__).resolve().parent.parent)
man = pd.read_csv(P.derived / "clone_manifest.csv")
feat = et.state_features(man, "BL")
Xc = [f"state_{s}" for s in et.STATES]

def design(unit, min_bl, extinct):
    g = et.growth_table(man, "BL", "EE", extinct=extinct)
    g = g[g.n_cells_src >= min_bl]
    if unit == "family_avg":
        g = g.groupby(["clone", "family"])["log2_growth"].mean().rename("log2_growth").reset_index()
    d = g.merge(feat.drop(columns="n_cells"), on=["clone", "family"]).dropna(subset=["log2_growth"])
    return d[Xc].to_numpy(), d["log2_growth"].to_numpy(), d["clone"].to_numpy()

rows = []
for unit, mb, ex in itertools.product(["per_passage", "family_avg"], [1, 5, 10, 20],
                                      ["pseudo", "censor"]):
    X, y, g = design(unit, mb, ex)
    s = et.grouped_cv_r(X, y, g, n_repeats=20, seed=0)
    rows.append(dict(unit=unit, min_bl_cells=mb, extinct=ex, n_rows=len(y),
                     n_clones=len(np.unique(g)), r_mean=s.mean(), r_sd=s.std()))
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(P.derived / "rung1_reconcile.csv", index=False)
