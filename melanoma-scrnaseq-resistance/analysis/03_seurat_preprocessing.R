# =============================================================================
# 03_seurat_preprocessing.R
#
# Builds the Seurat object from the QC-filtered counts, normalises, selects
# variable genes, scales them and runs PCA.
#
# Normalisation: LogNormalize (counts per 10,000, then log1p). GSE115978 is
# Smart-seq2 without UMIs, so library-size normalisation is an approximation;
# the published TPM matrix would be an alternative input.
# =============================================================================

library(Seurat)
library(dplyr)
library(readr)

set.seed(42)

proc_dir <- here::here("melanoma-scrnaseq-resistance", "data", "processed")

counts_filtered <- readRDS(file.path(proc_dir, "counts_filtered.rds"))
cell_annot <- read_csv(file.path(proc_dir, "cell_annotations.csv"), show_col_types = FALSE)

# The GEO annotation file names the cell id column "cells"
meta <- cell_annot %>%
  filter(cells %in% colnames(counts_filtered)) %>%
  as.data.frame()
rownames(meta) <- meta$cells
meta <- meta[colnames(counts_filtered), , drop = FALSE]

seu <- CreateSeuratObject(
  counts = counts_filtered,
  meta.data = meta,
  project = "GSE115978_melanoma_resistance",
  min.cells = 3,
  min.features = 200
)

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
# Scaling only the variable genes is enough for PCA and keeps memory low
seu <- ScaleData(seu, features = VariableFeatures(seu))
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 50, seed.use = 42)

saveRDS(seu, file.path(proc_dir, "seurat_preprocessed.rds"))
message(sprintf("Preprocessing done: %d genes x %d cells", nrow(seu), ncol(seu)))
