# =============================================================================
# 02_qc.R
#
# Cell quality control before building the Seurat object.
#
# GSE115978 is Smart-seq2 (full-length reads, no UMIs). QC uses the number of
# detected genes, total counts and mitochondrial percentage per cell.
# =============================================================================

library(dplyr)
library(readr)
library(Matrix)
library(ggplot2)

set.seed(42)

raw_dir  <- here::here("melanoma-scrnaseq-resistance", "data", "raw")
proc_dir <- here::here("melanoma-scrnaseq-resistance", "data", "processed")
fig_dir  <- here::here("melanoma-scrnaseq-resistance", "figures")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

counts <- read_csv(file.path(raw_dir, "GSE115978_counts.csv"), show_col_types = FALSE)
gene_names <- counts[[1]]
counts_mat <- as.matrix(counts[, -1])
rownames(counts_mat) <- gene_names
counts_mat <- Matrix(counts_mat, sparse = TRUE)

cell_annot <- read_csv(file.path(raw_dir, "GSE115978_cell.annotations.csv"), show_col_types = FALSE)

qc_metrics <- tibble(
  cell_id      = colnames(counts_mat),
  n_genes      = Matrix::colSums(counts_mat > 0),
  total_counts = Matrix::colSums(counts_mat)
)

mito_genes <- grep("^MT-", gene_names, value = TRUE)
qc_metrics$pct_mito <- if (length(mito_genes) > 0) {
  Matrix::colSums(counts_mat[mito_genes, , drop = FALSE]) / qc_metrics$total_counts * 100
} else {
  NA_real_
}

# Thresholds chosen from the distributions in the QC figures
MIN_GENES    <- 500
MIN_COUNTS   <- 1000
MAX_PCT_MITO <- 20

qc_metrics <- qc_metrics %>%
  mutate(pass_qc = n_genes >= MIN_GENES &
           total_counts >= MIN_COUNTS &
           (is.na(pct_mito) | pct_mito <= MAX_PCT_MITO))

message(sprintf("Cells: %d | passed QC: %d (%.1f%%)",
                nrow(qc_metrics), sum(qc_metrics$pass_qc), 100 * mean(qc_metrics$pass_qc)))

p_genes <- ggplot(qc_metrics, aes(x = n_genes, fill = pass_qc)) +
  geom_histogram(bins = 60) +
  geom_vline(xintercept = MIN_GENES, linetype = "dashed") +
  scale_fill_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "#D55E00"), name = "Passed QC") +
  labs(title = "Detected genes per cell", x = "Genes with counts > 0", y = "Cells") +
  theme_minimal()
ggsave(file.path(fig_dir, "qc_n_genes.png"), p_genes, width = 6, height = 4, bg = "white")

p_mito <- ggplot(qc_metrics, aes(x = pct_mito, fill = pass_qc)) +
  geom_histogram(bins = 60) +
  geom_vline(xintercept = MAX_PCT_MITO, linetype = "dashed") +
  scale_fill_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "#D55E00"), name = "Passed QC") +
  labs(title = "Mitochondrial counts per cell", x = "Mitochondrial counts (%)", y = "Cells") +
  theme_minimal()
ggsave(file.path(fig_dir, "qc_pct_mito.png"), p_mito, width = 6, height = 4, bg = "white")

cells_keep <- qc_metrics %>% filter(pass_qc) %>% pull(cell_id)
saveRDS(counts_mat[, cells_keep], file.path(proc_dir, "counts_filtered.rds"))
write_csv(qc_metrics, file.path(proc_dir, "qc_metrics.csv"))
write_csv(cell_annot, file.path(proc_dir, "cell_annotations.csv"))

message("QC done: ", proc_dir)
