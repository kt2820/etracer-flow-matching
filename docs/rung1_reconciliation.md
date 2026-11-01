# Reconciling the two rung-1 numbers (2026-11-01)

Rung 1 is the five baseline state fractions predicting BL→EE clone log2
growth on held-out clones (grouped 5-fold CV by clone, 20 repeats, ≥20 BL
cells per clone-family unit). Two sessions reported it differently:

| Session | r |
| --- | --- |
| 2026-10-30 | 0.323 ± 0.026 |
| 2026-10-31 (`nb01`, `nb02`) | 0.263 ± 0.025 |

The gap is about 2.4 sd, so it is not CV noise.

## What differed

1. **The clone table.** The Oct 30 run read `sc.metadata.BLEELE.csv` from
   Zenodo **16406856** (2025). This folder holds the newer deposit, Zenodo
   **19230230** (2026-03-26). Both tables cover the same 91,411 cells and the
   same state labels. Their clones map one-to-one with purity 1.000; only the
   IDs differ. The older table calls far more cells after baseline:

   | Stage | 19230230 (this folder) | 16406856 (Oct 30) |
   | --- | --- | --- |
   | BL | 47.5% | 48.5% |
   | EE | 47.0% | 90.7% |
   | LE | 56.6% | 95.4% |

2. **The growth definition.**
   - Oct 30 pooled the two EE passages into one fraction and floored
     extinction at a fraction of 1e-5.
   - Oct 31 averaged per-passage log ratios and floored extinction at 0.5
     cells.

## Result

`scripts/reconcile_clone_tables.py` crosses the two clone tables with every
growth definition (`derived/rung1_decomposition.csv`):

| Target | 19230230 | 16406856 | Share of target variance that is survive/die (19230230) |
| --- | --- | --- | --- |
| averaged, 0.5-cell floor (Oct 31) | 0.263 ± 0.025 | 0.267 ± 0.029 | 0.35 |
| pooled, 0.5-cell floor | 0.298 ± 0.020 | 0.296 ± 0.022 | 0.52 |
| averaged, 1e-5 floor | 0.349 ± 0.017 | 0.339 ± 0.019 | 0.66 |
| pooled, 1e-5 floor (Oct 30) | 0.344 ± 0.014 | 0.329 ± 0.023 | 0.81 |
| **survival alone (binary)** | **0.330 ± 0.020** | **0.331 ± 0.020** | 1.00 |
| growth among survivors only | 0.145 ± 0.033 | 0.140 ± 0.025 | — |

The Oct 30 recipe on the old table gives 0.329. The remaining 0.006 to 0.323
comes from inner-CV details: MSE vs Pearson for penalty selection, and the
penalty grid.

## Conclusions

- **The clone table barely matters for rung 1.** The two tables agree within
  0.015 under every definition, because baseline calls are almost identical.
  The older table's extra EE/LE calls would matter for anything measured after
  baseline, such as descendant expression.
- **The growth definition sets the number.** The extinction floor is an
  arbitrary constant, and it decides how far extinct clones sit below
  survivors. A deeper floor makes the target more and more a survive/die
  indicator, and r rises toward the survival-only score. Neither 0.323 nor
  0.263 is wrong; they measure different mixtures of survival and growth.
- **The predictable signal is survival.** Survival alone scores 0.330.
  Growth among survivors scores 0.145. This backs the reformulation: clone
  survival as the primary (binary, calibrated) target, and log growth among
  survivors as the secondary target.
- **Earlier expression rungs need re-scoring.** The Oct 30 pseudobulk (0.376)
  and kernel mean embedding (0.441) used the Oct 30 target and table. Re-score
  them under the chosen target before comparing them with anything in
  `nb01`/`nb02`. `derived/headroom.csv` currently mixes the two setups.
