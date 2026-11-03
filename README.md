# Lineage-coupled flow matching on eTRACER

My course project on optimal transport and flow matching.

**Problem statement.** Using the eTRACER dataset I'm trying to understand
how unbalanced flow-matching models of single-cell tumor evolution behave
when cells are coupled by clonal barcodes vs an inferred OT coupling. I will
explore how the different couplings shape clone growth and cell expression
under immune pressure. And I will explore how different ways of modeling
growth affect the learned dynamics.

**The data.** eTRACER (Yang, Hou, Wang et al., bioRxiv 2025.08.06.668639;
Zenodo record 19230230): barcoded mouse lung-cancer clones grown as tumors
under CD8 T-cell pressure and passaged three generations deep, with
single-cell RNA-seq at every stage.

## Notebooks (`data_exploration/`)

These are meant to be read in the following order:

1. **`etracer_overview.ipynb`.** A quick overview of the dataset and conventions.
2. **`survival_motivation.ipynb`.** Trying to understand how a
   clone's Primary-tumor cells could predict its survival.
3. **`coupling_ablation.ipynb`.** Motivation for the proposed experiment.

`GLOSSARY.md` defines the terms used across the notebooks and points to
where each is computed.

The notebooks import `src/etracer.py` and read cached arrays from `derived/`
and `review/`, which are not committed for now.

## Data

The count matrix the notebooks start from, `derived/etracer_counts.h5ad`
(91,411 cells x 32,285 genes, raw integer UMI counts converted from the
deposit's Seurat object), is available on
[Google Drive](https://drive.google.com/file/d/1ptOXFz8a9KazCyosHFQ5SWC4ccRepFjs/view?usp=drive_link).
Download it and place it at `derived/etracer_counts.h5ad`.
