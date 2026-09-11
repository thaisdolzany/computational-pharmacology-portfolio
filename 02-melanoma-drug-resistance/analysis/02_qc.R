# =============================================================================
# 02_qc.R
#
# Controle de qualidade das células antes da criação do objeto Seurat.
#
# Nota: GSE115978 é Smart-seq2 (full-length), não 10x/UMI. QC aqui é
# baseado em número de genes detectados, contagem total e % mitocondrial
# por célula.
# =============================================================================

library(dplyr)
library(readr)
library(Matrix)
library(ggplot2)

set.seed(42)

raw_dir  <- here::here("02-melanoma-drug-resistance", "data", "raw")
proc_dir <- here::here("02-melanoma-drug-resistance", "data", "processed")
fig_dir  <- here::here("02-melanoma-drug-resistance", "figures")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Carrega matriz de contagens e anotações de célula
# -----------------------------------------------------------------------------
counts <- read_csv(file.path(raw_dir, "GSE115978_counts.csv"))
gene_col <- names(counts)[1]
gene_names <- counts[[gene_col]]
counts_mat <- as.matrix(counts[, -1])
rownames(counts_mat) <- gene_names
counts_mat <- Matrix(counts_mat, sparse = TRUE)

cell_annot <- read_csv(file.path(raw_dir, "GSE115978_cell.annotations.csv"))

# -----------------------------------------------------------------------------
# Métricas de QC por célula
# -----------------------------------------------------------------------------
qc_metrics <- tibble(
  cell_id       = colnames(counts_mat),
  n_genes       = Matrix::colSums(counts_mat > 0),
  total_counts  = Matrix::colSums(counts_mat)
)

mito_genes <- grep("^MT-", gene_names, value = TRUE)
if (length(mito_genes) > 0) {
  qc_metrics$pct_mito <- Matrix::colSums(counts_mat[mito_genes, , drop = FALSE]) /
    qc_metrics$total_counts * 100
} else {
  qc_metrics$pct_mito <- NA_real_
}

# -----------------------------------------------------------------------------
# Filtros de QC (thresholds iniciais — revisar após inspeção visual)
# -----------------------------------------------------------------------------
MIN_GENES <- 500
MIN_COUNTS <- 1000
MAX_PCT_MITO <- 20

qc_metrics <- qc_metrics %>%
  mutate(
    pass_qc = n_genes >= MIN_GENES &
      total_counts >= MIN_COUNTS &
      (is.na(pct_mito) | pct_mito <= MAX_PCT_MITO)
  )

message(sprintf(
  "Células totais: %d | Passaram QC: %d (%.1f%%)",
  nrow(qc_metrics), sum(qc_metrics$pass_qc),
  100 * mean(qc_metrics$pass_qc)
))

# -----------------------------------------------------------------------------
# Plots de QC
# -----------------------------------------------------------------------------
p_genes <- ggplot(qc_metrics, aes(x = n_genes, fill = pass_qc)) +
  geom_histogram(bins = 60) +
  geom_vline(xintercept = MIN_GENES, linetype = "dashed") +
  labs(title = "Genes detectados por célula", x = "n_genes", y = "count") +
  theme_minimal()
ggsave(file.path(fig_dir, "qc_n_genes.png"), p_genes, width = 6, height = 4)

p_mito <- ggplot(qc_metrics, aes(x = pct_mito, fill = pass_qc)) +
  geom_histogram(bins = 60) +
  geom_vline(xintercept = MAX_PCT_MITO, linetype = "dashed") +
  labs(title = "% mitocondrial por célula", x = "% mito", y = "count") +
  theme_minimal()
ggsave(file.path(fig_dir, "qc_pct_mito.png"), p_mito, width = 6, height = 4)

# -----------------------------------------------------------------------------
# Salva matriz filtrada + metadados para o próximo passo
# -----------------------------------------------------------------------------
cells_keep <- qc_metrics %>% filter(pass_qc) %>% pull(cell_id)
counts_filtered <- counts_mat[, cells_keep]

saveRDS(counts_filtered, file.path(proc_dir, "counts_filtered.rds"))
write_csv(qc_metrics, file.path(proc_dir, "qc_metrics.csv"))
write_csv(cell_annot, file.path(proc_dir, "cell_annotations.csv"))

message("QC concluído. Matriz filtrada salva em: ", proc_dir)
