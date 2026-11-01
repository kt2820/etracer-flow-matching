#!/usr/bin/env Rscript
# Convert scRNA/sc_UMI_count.rds to a portable sparse matrix.
#
# RUN THIS ON GREENE, not through the Cowork VM: that shell has no R, and the
# local disk was ~97% full (18 GB free) when this was written. The .rds is
# 1.0 GB compressed and expands to roughly 5 GB in memory.
#
#   module load r/4.3.1
#   Rscript convert_seurat_to_h5ad.R <in.rds> <outdir>
#
# Writes MatrixMarket + two text files, which `assemble_h5ad.py` turns into an
# h5ad. MatrixMarket is chosen over SaveH5Seurat on purpose: no SeuratDisk
# dependency, and it works whether the .rds holds a Seurat object or a bare
# dgCMatrix.

args <- commandArgs(trailingOnly = TRUE)
in_rds <- if (length(args) >= 1) args[1] else "sc_UMI_count.rds"
outdir <- if (length(args) >= 2) args[2] else "derived/matrix"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(Matrix))

message("reading ", in_rds, " ...")
obj <- readRDS(in_rds)
message("class: ", paste(class(obj), collapse = ", "))

counts <- if (inherits(obj, "Seurat")) {
  suppressPackageStartupMessages(library(Seurat))
  # Seurat v5 keeps counts in a layer; v4 in slot="counts".
  tryCatch(
    SeuratObject::LayerData(obj, assay = "RNA", layer = "counts"),
    error = function(e) Seurat::GetAssayData(obj, assay = "RNA", slot = "counts")
  )
} else if (inherits(obj, "dgCMatrix") || inherits(obj, "Matrix")) {
  obj
} else if (inherits(obj, "matrix")) {
  as(obj, "dgCMatrix")
} else {
  stop("unrecognised object class: ", paste(class(obj), collapse = ", "))
}

message("dim: ", nrow(counts), " genes x ", ncol(counts), " cells")
message("nnz: ", length(counts@x))
stopifnot(all(counts@x == floor(counts@x)))   # must still be raw integer UMI

writeMM(counts, file.path(outdir, "counts.mtx"))
writeLines(rownames(counts), file.path(outdir, "genes.txt"))
writeLines(colnames(counts), file.path(outdir, "cells.txt"))

# Per-cell totals, so the Python side can re-verify against nCount_RNA.
tot <- Matrix::colSums(counts)
write.csv(data.frame(cell = colnames(counts), total = tot),
          file.path(outdir, "cell_totals.csv"), row.names = FALSE)

message("done -> ", outdir)
message("NOTE: counts.mtx is plain text and will be several GB. gzip it ",
        "before copying, or run the expression notebooks on Greene directly.")
