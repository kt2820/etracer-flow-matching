"""Chart the conditioning ladder: docs/ladder.png from derived/ladder_*.csv."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
S = pd.read_csv(ROOT / "derived/ladder_scores.csv")
C = pd.read_csv(ROOT / "derived/ladder_ceiling.csv").set_index("target")
Pq = pd.read_csv(ROOT / "derived/ladder_paired.csv").set_index("target")

SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
RAMP = ["#86b6ef", "#2a78d6", "#104281"]          # ordinal blue: richer description -> darker
FEATS = ["state fractions (5)", "pseudobulk mean (400 genes)", "kernel mean embedding"]
LABEL = ["State\nfractions", "Pseudobulk\nmean", "Distribution\n(kernel embedding)"]
PANELS = [("survival (primary)", "Primary target: does the clone survive to early evasion?"),
          ("growth | survived", "Secondary target: log growth among survivors")]

plt.rcParams.update({"font.size": 10, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK})
fig, axes = plt.subplots(1, 2, figsize=(11, 5.2), sharey=True, facecolor=SURF)
for ax, (tgt, title) in zip(axes, PANELS):
    ax.set_facecolor(SURF)
    d = S[S.target == tgt].set_index("features").loc[FEATS]
    c = C.loc[tgt]
    ax.axhspan(c.ceiling_lo, c.ceiling_hi, color=GRID, alpha=0.35, lw=0, zorder=0)
    ax.axhline(c.ceiling, color=INK2, lw=1.2, ls="--", zorder=1)
    ax.text(2.45, c.ceiling + 0.012, f"reliability ceiling {c.ceiling:.2f}", ha="right", va="bottom",
            fontsize=8.5, color=INK2)
    nm = d.null_mean.mean(); ns = d.null_sd.mean()
    ax.axhspan(nm - ns, nm + ns, color="#d6d5cf", alpha=0.5, lw=0, zorder=0)
    x = np.arange(3)
    ax.bar(x, d.r_mean, width=0.58, color=RAMP, edgecolor=SURF, linewidth=2, zorder=2)
    ax.errorbar(x, d.r_mean, yerr=d.r_sd, fmt="none", ecolor=INK, elinewidth=1.2, capsize=4, zorder=3)
    for xi, (m, s) in enumerate(zip(d.r_mean, d.r_sd)):
        ax.text(xi, m + s + 0.012, f"{m:.2f}", ha="center", va="bottom", fontsize=10, color=INK, weight="bold")
    ax.set_xticks(x, LABEL); ax.tick_params(axis="x", length=0)
    ax.set_title(title, loc="left", fontsize=10.5, color=INK, pad=10)
    ax.grid(axis="y", color=GRID, lw=0.8, zorder=0); ax.set_axisbelow(True)
    for sp_ in ["top", "right", "left"]: ax.spines[sp_].set_visible(False)
    n = int(d.n_units.iloc[0]); k = int(d.n_clones.iloc[0])
    note = f"{n} clone-family units, {k} clones  ·  ceiling 95% CI {c.ceiling_lo:.2f}–{c.ceiling_hi:.2f}"
    if tgt.startswith("survival"):
        p = Pq.loc[tgt]
        note += (f"  ·  AUC {d.auc_mean.iloc[0]:.2f} / {d.auc_mean.iloc[1]:.2f} / {d.auc_mean.iloc[2]:.2f}"
                 f"\ndistribution beats mean in {int(p.wins)} of {int(p.of)} repeats (+{p.kme_minus_mean:.3f})")
    ax.text(0, -0.30, note, transform=ax.transAxes, fontsize=8.5, color=INK2, va="top")
axes[0].set_ylabel("Pearson r on held-out clones")
axes[0].set_ylim(-0.12, 0.85)
fig.suptitle("Baseline expression predicts which clones survive immune pressure, not how fast survivors grow",
             x=0.01, ha="left", fontsize=12, weight="bold", color=INK)
fig.text(0.01, 0.935, "eTRACER, baseline → early evasion · ridge on baseline-cell features · grouped 5-fold CV by clone, mean ± sd of 20 repeats\n"
         "Dashed line: reliability ceiling (same clone, other founder families), shaded = 95% CI · grey band at zero: permutation null ±1 sd",
         fontsize=9, color=INK2, va="top")
fig.tight_layout(rect=(0, 0.02, 1, 0.86))
fig.savefig(ROOT / "docs/ladder.png", dpi=200, facecolor=SURF)
print("wrote docs/ladder.png")
