# Lineage-coupled flow matching on eTRACER

Course project for a graduate class on optimal transport and flow matching.

**Problem statement.** Using the eTRACER dataset I'm trying to understand
how unbalanced flow-matching models of single-cell tumor evolution behave
when cells are coupled by clonal barcodes vs an inferred OT coupling. I will
explore how the inferred coupling shapes clone growth and cell expression
under immune pressure. And I will explore how different ways of modeling
growth affect the learned dynamics.

**The data.** eTRACER (Yang, Hou, Wang et al., bioRxiv 2025.08.06.668639;
Zenodo record 19230230): barcoded mouse lung-cancer clones grown as tumors
under CD8 T-cell pressure and passaged three generations deep, with
single-cell RNA-seq at every stage. Flow matching between two tumors needs a
coupling: which Primary cells are paired with which Secondary cells. Usually
that coupling is inferred by optimal transport on expression. Here the
barcodes give it exactly, clone by clone, including which clones grew and
which were lost.

## Notebooks (`data_exploration/`)

Read them in this order. Each one ends where the next one begins.

1. **`etracer_overview.ipynb` — what the data is.** The experimental
   design (3 founder families, 3 → 6 → 6 tumors), how the deposit's names
   map to ours, the five tumor states, the 173 clones and their trees, and
   the clonal bottleneck: most clones are never seen again after the
   Primary tumor, and survival rises with clone size. Runs on the metadata
   and trees alone.
2. **`survival_motivation.ipynb` — selection acts on clones.** Does a
   clone's Primary-tumor cells predict its survival beyond its size? Fate
   repeats across independent mice once size is accounted for, and
   survivors already differ in state mix and expression before selection,
   though mostly through the state mix. So which cells have descendants is
   a clone-level fact that expression reflects only in part. Needs the
   count matrix.
3. **`coupling_ablation.ipynb` — what the barcodes buy.** Flow matching
   needs a coupling between Primary and Secondary cells. Compares the exact
   lineage coupling (cells paired within clone, clone mass from the
   observed growth) with the coupling standard unbalanced OT would infer
   from expression alone, on ancestry, growth, state transitions, velocity
   targets and generalisation to held-out clones. Ends with the ranked next
   steps for training the flows.

`GLOSSARY.md` defines the terms used across the notebooks and points to
where each is computed.

The notebooks import `src/etracer.py` and read cached arrays from `derived/`
and `review/`, which are not committed.

## Data

The count matrix the notebooks start from, `derived/etracer_counts.h5ad`
(91,411 cells x 32,285 genes, raw integer UMI counts converted from the
deposit's Seurat object), is available on
[Google Drive](https://drive.google.com/file/d/1ptOXFz8a9KazCyosHFQ5SWC4ccRepFjs/view?usp=drive_link).
Download it and place it at `derived/etracer_counts.h5ad`.
