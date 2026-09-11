# =============================================================================
# 04_seurat_preprocessing.R
#
# Cria o objeto Seurat a partir da matriz já log2-normalizada (log2[(TPM/10)+1]),
# identifica genes variáveis, escala e roda PCA.
#
# Nota: como os dados já chegam log-transformados (não são contagens
# brutas), NÃO reaplicamos NormalizeData/LogNormalize (que assumiria
# contagens tipo-UMI). Em vez disso, copiamos a matriz para o slot "data"
# diretamente.
# =============================================================================

library(Seurat)
library(dplyr)
library(readr)

set.seed(42)

proc_dir <- here::here("01-melanoma-single-cell-biomarkers", "data", "processed")

log_tpm_filtered <- readRDS(file.path(proc_dir, "log_tpm_filtered.rds"))
qc_metrics <- read_csv(file.path(proc_dir, "qc_metrics.csv"))

# -----------------------------------------------------------------------------
# Monta metadados alinhados às células filtradas
# -----------------------------------------------------------------------------
meta <- qc_metrics %>%
  filter(pass_qc, cell_id %in% colnames(log_tpm_filtered)) %>%
  as.data.frame()
rownames(meta) <- meta$cell_id
meta <- meta[colnames(log_tpm_filtered), , drop = FALSE]

# -----------------------------------------------------------------------------
# Cria o objeto Seurat (usando a matriz log-transformada como "counts" só
# para satisfazer o construtor; corrigido logo abaixo)
# -----------------------------------------------------------------------------
seu <- CreateSeuratObject(
  counts = log_tpm_filtered,
  meta.data = meta,
  project = "GSE120575_melanoma_biomarkers",
  min.cells = 3,
  min.features = 200
)

# -----------------------------------------------------------------------------
# Dados já normalizados/log-transformados na origem: copia para o slot
# "data" em vez de rodar NormalizeData (que assumiria contagens brutas)
# -----------------------------------------------------------------------------
seu <- SetAssayData(seu, slot = "data", new.data = GetAssayData(seu, slot = "counts"))

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
# Salva objeto para o próximo passo (clustering + biomarcadores)
# -----------------------------------------------------------------------------
saveRDS(seu, file.path(proc_dir, "seurat_preprocessed.rds"))

message("Pré-processamento concluído. Objeto Seurat salvo em: ", proc_dir)
message(sprintf("Dimensões finais: %d genes x %d células", nrow(seu), ncol(seu)))
