# Glossary for the data_exploration notebooks

Definitions of the non-trivial terms and ideas used in `etracer_overview`,
`survival_motivation` and `coupling_ablation`. Each entry says where the
thing is defined or computed in the code, so a term can be read next to the
lines that use it. Notebook code comments point back here as
`GLOSSARY.md: term`.

Notebooks are abbreviated **OV** (`etracer_overview.ipynb`),
**SM** (`survival_motivation.ipynb`) and **CA** (`coupling_ablation.ipynb`).

---

## 1. The experiment and its names

**eTRACER.** The dataset: mouse lung-cancer cells (EGFR-L858R; Trp53-null)
carrying lineage barcodes, grown as tumors under CD8 T-cell pressure and
passaged three generations deep (Yang, Hou, Wang et al., bioRxiv
2025.08.06.668639). Zenodo record 19230230.
*Code:* OV, opening markdown and the "How the deposit's names map to ours" cell.

**OT-I T cells.** CD8 T cells engineered to recognise one antigen (ovalbumin
peptide) that the tumor cells express. Transferring them into the mouse
applies immune selection on the tumor. Primary tumors grew without them;
Secondary and Tertiary tumors grew under OT-I transfer.
*Code:* OV, opening markdown.

**Family (A / B / C).** One founder (Primary) tumor and everything seeded
from it. The deposit calls them a3026 / a3028 / a3030. The three families
were seeded from the same barcoded pool, so the same clone can appear in
several families.
*Code:* `FAMILY_OF` in OV, SM and CA.

**Generation / stage (Primary / Secondary / Tertiary).** How many
transplants a tumor is from the founder. The deposit and paper call them
BL / EE / LE (baseline, early evasion, late evasion).
*Code:* `STAGE_OF` in OV and SM; `TUMOR_OF` in CA (which also splits by branch).

**Tumor name (A1, A2a, A3b, …).** Letter = family, digit = generation,
a/b = which of the two Secondary branches. A3a grows from A2a. CA writes
the two Secondary tumors of a family as "Secondary a" and "Secondary b".
*Code:* `tumor_name()` in OV and SM; `parent_of()` in OV; `TUMOR_OF` in CA.

**Branch.** The deposit's passage index inside a family: 0 for Primary,
1 / 2 for the two Secondary tumors, 1_1 / 2_1 for the Tertiary tumor grown
from each. In CA it is the key that tells Secondary a from Secondary b.
*Code:* CA, `TUMOR_OF` cell.

**Sister tumors.** The two Secondary tumors of one family: separate mice,
each engrafted from half of the same Primary tumor. Because they are
independent engraftments, one can serve as the measurement the other is
scored against.
*Code:* CA §2, `SISTER`; CA §3 reference rows.

**Siblings, not ancestors.** Half of each tumor was transplanted and the
other half sequenced, so the sequenced cells are siblings of the cells that
seeded the next tumor, never their direct ancestors. A coupling between
Primary and Secondary cells is therefore a clone-to-clone statement, not
cell-level ancestry.
*Code:* OV opening markdown; CA "Couplings" markdown.

---

## 2. Lineage

**Static barcode.** A fixed DNA barcode inserted before the experiment.
Cells sharing a static barcode descend from the same founder cell.
*Code:* OV §3 markdown.

**Clone.** The set of cells sharing static barcodes; numbered "clone 1",
"clone 67", … so they never clash with tumor names. 173 clones in total.
*Code:* `clone_name()` in OV; `meta["clone"]` everywhere.

**Clone call / called.** Whether a cell's barcode was read well enough to
assign it to a clone. About half of cells are called; the rest have a
transcriptome but no lineage label.
*Code:* `meta["called"]`; `call_rate` table in OV §2.

**Unit (clone-in-family, clone-family pair).** One clone in one family's
tumor — the row on which every statistic is computed. Clone 12 in family A
and clone 12 in family C are two units.
*Code:* `units` in SM ("The unit" section); the per-dataset clone tables in CA.

**Evolving barcode / Cas9 edits.** A second barcode that accumulates Cas9
edits as cells divide, so cells that are more closely related share more
edits. Used to build trees inside each clone.
*Code:* OV §4 markdown.

**Clone tree (phylogeny).** An IQ-TREE tree over one clone's cells, built
from the evolving edits, with cells from all generations and families in the
same tree. Leaves are cells. Branch lengths are edit distance, not time.
*Code:* OV §4, `trees`, `tree_summary`, `plot_clone_tree()`.

**Newick.** The text format the trees are stored in. `et.read_newick`
parses it; `et.tree_layout` assigns x/y positions for drawing.
*Code:* OV §4.

---

## 3. Expression

**UMI (unique molecular identifier).** A tag on each captured RNA molecule.
The UMI count of a cell is how many molecules were sequenced from it, a
proxy for sequencing depth and cell size.
*Code:* `n_umi` in SM.

**Sequencing depth.** The number of molecules read per cell. Deeper cells
show more non-zero genes after normalisation, which can mimic biology, so it
is adjusted for.
*Code:* SM, `depth` and `Z_depth` ("Sequencing depth" cell).

**log-CPM / log1p(CP10k).** Counts per 10,000 molecules, then log(1 + x).
Puts cells of different depth on a comparable scale.
*Code:* CA, `Lc` array loaded from `flow_v1_cells.npz`.

**Pseudobulk.** The mean (or sum) of expression over a group of cells,
treated as one sample. Here: the mean log-normalised expression of a unit's
Primary cells.
*Code:* `et.pseudobulk()` call in SM ("Pseudobulk" cell).

**Tumor state (the five states).** Clusters the eTRACER preprint found and
named: Hypoxic, Proliferative, EMT, Metabolic, Interferon response. Each
cell carries one. The deposit numbers them 1–5.
*Code:* `STATE_NAMES`, `STATE_OF` in OV and SM; `state_name` in CA.

**State mix / state fractions.** For a unit, the share of its cells in each
of the five states. Five numbers that sum to one.
*Code:* `mix` in SM §2; `comp` in OV §2.

**Gene panel.** A fixed set of genes chosen without looking at clone fate,
so a model cannot be helped by gene selection that peeked at the answer.
SM uses genes expressed in at least 5% of Primary cells; CA uses the 2,371
nuclear (non-mitochondrial) genes of the cached panel.
*Code:* SM `genes`; CA `nuclear` mask.

**PCs (principal components).** A 50-dimensional linear compression of the
gene panel in which the couplings are computed, fitted on the Primary and
Secondary cells used.
*Code:* CA, PCA cell (`V`, `Z`).

**Harmony UMAP.** The deposit's 2-D embedding of all cells, after Harmony
batch correction. Used only for pictures.
*Code:* OV §2, `umap1` / `umap2`.

---

## 4. Fate and selection

**Survival (of a clone).** A clone in a parent tumor survived if at least one
of its cells is seen in the child tumor. In SM the two Secondary tumors are
pooled: seen in either counts.
*Code:* OV §5 `surv`; SM `units["survived"]`.

**Lost / extinct vs missed.** A clone not seen in the child is "lost". It may
be extinct or merely undetected, because only about half of cells are
called. Small clones are missed more often.
*Code:* OV §5 "Survival by how big the clone was" cell.

**Clone growth.** A clone's share of the child tumor relative to its share
of the parent, on a log2 scale. In CA, the growth a coupling *implies* for a
clone is its source marginal (mass leaving the clone's Primary cells)
divided by the uniform weight, averaged over the clone's cells.
*Code:* CA §2, `g_here`, `g_sister`, `g_inf`.

**Bottleneck.** The sharp fall in clonal diversity at each transplant: a few
clones take over.
*Code:* OV §3 "Clone size distribution" and "Effective number of clones" cells.

**Effective number of clones.** exp(Shannon entropy) of clone shares: the
number of equally sized clones that would give the same diversity.
*Code:* OV §3 `effective_clones()`.

**Beating the size prediction (`beat_size`).** Observed survival (0/1) minus
the probability of survival predicted from clone size and family. Positive:
the clone did better than clones of its size usually do.
*Code:* SM §1, `beat_size`; `fit_logistic()`.

**Fate is a property of the clone.** The claim that, after accounting for
size, a clone that survives more than expected in one family does so in
others. Tested by correlating `beat_size` for the same clone across families
against a within-family shuffle.
*Code:* SM `cross_family_r()`, `shuffle_within_family()`.

**Selection read as state change.** What a lineage-free coupling does when
a state's cells mostly die: the mass has to go somewhere, so it is moved to
other states and looks like a transition. CA's Metabolic row is the example.
*Code:* CA §3, the Metabolic / Interferon transition rows.

---

## 5. Couplings and the model

**Flow matching.** A way to learn a vector field $v(x, t)$ that moves
samples from one distribution to another by following $dx/dt = v$. Training
regresses $v$ onto the straight-line velocity between paired source and
target points, so the pairing — the coupling — decides what is learned.
*Code:* CA "Couplings" markdown.

**Coupling.** A joint distribution $\pi$ over (Primary cell, Secondary cell)
pairs. Its row sums say how much mass leaves each Primary cell; its columns
say where it lands. Everything CA measures is read off $\pi$.
*Code:* CA `coupling_summaries()`.

**Lineage coupling.** Cells are paired only within their clone (static
barcodes, exact). Inside a clone the pairing is balanced entropic OT; the
clone's block carries its observed share of the Secondary tumor, so growth
and loss are the measured ones. Clones absent at Secondary carry zero mass.
*Code:* CA, "Lineage" block of the couplings loop; `out["lineage"]`.

**Inferred coupling.** Standard unbalanced entropic OT between all Primary
and all Secondary cells, with no lineage: uniform marginals, a KL penalty on
the source side so mass can grow or die. The UOT-FM default (Eyring et al.
2024).
*Code:* CA, "Inferred" block; `out["inferred τ=…"]`.

**Entropic OT / Sinkhorn.** OT with an entropy term of weight $\varepsilon$,
solved by alternately rescaling rows and columns of $K = \exp(-C/\varepsilon)$.
Cost is squared Euclidean distance in PC space divided by its median, so
$\varepsilon$ means the same thing in every dataset.
*Code:* CA `sinkhorn()`, `EPS`, `sqdist()`.

**Unbalanced OT (τ).** Lets the source marginal deviate from uniform at a
KL cost of weight $\tau$: small $\tau$ is nearly free growth, large $\tau$ is
nearly balanced. CA runs three values, since growth is the thing lineage
supplies and the inferred coupling has to guess it.
*Code:* CA `TAUS`; the `tau` argument of `sinkhorn()`.

**Ancestry (within-clone mass, top-1 clone).** Two readings of whether a
coupling respects clones: the share of transported mass that stays inside
the true clone, and the share of Secondary cells whose largest-mass source
clone is their own. Chance levels come from clone shares.
*Code:* CA §1, `ancestry` table.

**State transition matrix.** $T_{ij}$: the share of transported mass from
Primary state $i$ to Secondary state $j$, row-normalised so each row is
"where does state $i$ go". Rows of two couplings are compared by total
variation distance, against the variation between sister tumors as reference.
*Code:* CA §3, `TMATS`, `row_norm()`, `transitions`.

**Velocity target / displacement.** The flow-matching regression target for
a Primary cell: its expected position under the coupling minus where it
is, $\mathbb{E}_\pi[x_1 \mid x_0] - x_0$. CA compares the two couplings'
displacements per cell (cosine, length ratio) and splits the lineage
displacement into the clone's mean shift and within-clone redistribution.
*Code:* CA §4, `disp` in `coupling_summaries()`, `velocity` table.

**Carry-forward.** Predicting that a clone's Secondary mean equals its
Primary mean. The baseline a learned map has to beat; a constant tumor-wide
shift adds nothing after centring and so scores the same.
*Code:* CA §5, `carry-forward` and `tumor shift` columns.

---

## 6. Statistics and validation

**Grouped cross-validation by clone.** Folds are made so all cells or units
of a clone fall in the same fold. A clone is never in both training and
test, because the same clone in two families would otherwise leak.
*Code:* `et.grouped_folds()` in CA §5.

**Ridge regression.** Linear regression with an L2 penalty; the penalty is
chosen by an inner cross-validation on training clones. In CA it is the
one-step linear "flow" from a cell's position to its displacement.
*Code:* CA §5 `ridge_closed()`, `ridge_apply()`, `LAM`.

**Logistic regression / IRLS.** Fits P(survive) as a logistic function of
covariates. IRLS (iteratively reweighted least squares) is Newton's method
on the log-likelihood.
*Code:* SM `fit_logistic()`.

**AUC.** Area under the ROC curve: the probability that a random survivor is
scored above a random lost clone. 0.5 is chance. Computed as a rank statistic.
*Code:* SM `auc()`.

**Pearson r / Spearman.** Linear and rank correlation. Used across units
(prediction vs outcome), across clones (implied vs measured growth), and
pooled across (clone, PC) entries in CA §5.
*Code:* `et.pearson()`, `et.spearman()`.

**Total variation distance.** $\tfrac{1}{2}\sum_j |p_j - q_j|$ between two
rows of a transition matrix: 0 when identical, 1 when disjoint.
*Code:* CA §3, `tv`.

**Residualise / adjust for covariates.** Regress a variable on covariates
(size, family, depth, …) and keep the residual. Comparing residuals is a
comparison at matched covariates.
*Code:* SM `residualise()`, `Z_size`, `Z_depth`, `Z_state`.

**Permutation null.** Repeat the analysis with the outcome labels shuffled
(here: give each clone another clone's outcomes, all its families at once)
to see what the statistic does when there is nothing to find.
*Code:* SM `survival_gene_test()`; `shuffle_within_family()`.

**Pooled null across genes.** One null distribution built from all genes'
permuted statistics together, so the smallest attainable p is
1 / (permutations × genes).
*Code:* SM `survival_gene_test()`, `pooled`.

**FDR / Benjamini–Hochberg (q-value).** Controls the expected share of false
positives among the genes called significant. q < 0.10 means about 10% of
the called genes are expected to be false.
*Code:* SM `survival_gene_test()`, `q`.

**Bootstrap over clones.** Resample units with replacement and recompute the
statistic to get a confidence interval.
*Code:* SM `boot_ci()` and the state-mix bootstrap in §2.

**Ceiling.** How well a measurement reproduces itself between independent
mice; an upper bound on what any predictor of it could reach. In CA §2 the
lineage coupling's growth correlation with the sister tumor is that ceiling.
*Code:* CA §2, `lineage_r`.

**Confounded.** Two or more causes that always change together in the data,
so their effects cannot be separated. Here: T-cell pressure, passage, and
elapsed time all change at the Primary → Secondary step, and there is no
untreated arm.
*Code:* OV opening markdown; CA "Next steps" row 3.
