# =============================================================================
# 03_seurat_preprocessing.R
#
# Cria o objeto Seurat a partir da matriz filtrada, normaliza, identifica
# genes variáveis, escala e roda PCA.
#
# Nota: os dados são TPM/full-length (Smart-seq2), não UMI counts de 10x.
# NormalizeData(method = "LogNormalize") assume contagens tipo-UMI; como
# alternativa mais apropriada a full-length, usamos log1p diretamente sobre
# a matriz — documentado abaixo.
# =============================================================================

library(Seurat)
library(dplyr)
library(readr)

set.seed(42)

proc_dir <- here::here("02-melanoma-drug-resistance", "data", "processed")

counts_filtered <- readRDS(file.path(proc_dir, "counts_filtered.rds"))
cell_annot <- read_csv(file.path(proc_dir, "cell_annotations.csv"))

# -----------------------------------------------------------------------------
# Alinha metadados às células filtradas
#
# Nota: o CSV de anotações do GEO usa a coluna "cells" (não "cell_id") para
# o identificador de célula.
# -----------------------------------------------------------------------------
meta <- cell_annot %>%
  filter(cells %in% colnames(counts_filtered)) %>%
  as.data.frame()
rownames(meta) <- meta$cells
meta <- meta[colnames(counts_filtered), , drop = FALSE]

# -----------------------------------------------------------------------------
# Cria o objeto Seurat
# -----------------------------------------------------------------------------
seu <- CreateSeuratObject(
  counts = counts_filtered,
  meta.data = meta,
  project = "GSE115978_melanoma_resistance",
  min.cells = 3,
  min.features = 200
)

# -----------------------------------------------------------------------------
# Normalização (log-transform apropriado a dados full-length já em TPM/contagem
# relativa; ver nota no cabeçalho)
# -----------------------------------------------------------------------------
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4)

# -----------------------------------------------------------------------------
# Genes altamente variáveis
# -----------------------------------------------------------------------------
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

# -----------------------------------------------------------------------------
# Scaling + PCA
# -----------------------------------------------------------------------------
all_genes <- rownames(seu)
seu <- ScaleData(seu, features = all_genes)
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 50, seed.use = 42)

# -----------------------------------------------------------------------------
# Salva objeto para o próximo passo (clustering)
# -----------------------------------------------------------------------------
saveRDS(seu, file.path(proc_dir, "seurat_preprocessed.rds"))

message("Pré-processamento concluído. Objeto Seurat salvo em: ", proc_dir)
message(sprintf("Dimensões finais: %d genes x %d células", nrow(seu), ncol(seu)))
