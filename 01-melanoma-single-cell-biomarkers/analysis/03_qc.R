# =============================================================================
# 03_qc.R
#
# Carrega a matriz TPM (genes x células), junta com os metadados por
# célula extraídos em 02_parse_metadata.R, e aplica QC.
#
# Notas importantes sobre este dataset:
#  - O arquivo de matriz tem 2 linhas de cabeçalho: linha 1 = ID da célula,
#    linha 2 = ID do paciente/timepoint (redundante com cell_metadata.csv,
#    ignorada aqui).
#  - Os valores já vêm como log2[(TPM/10)+1] (não são contagens brutas).
#    Como log2(0/10 + 1) = 0, o valor 0 no arquivo corresponde exatamente a
#    TPM = 0 — a estrutura esparsa da matriz continua válida.
#  - Não há "no.of.reads"/"no.of.genes" pré-calculados pelos autores para
#    este dataset (diferente do GSE115978 usado no Projeto 02) — as
#    métricas de QC abaixo são calculadas a partir da própria matriz.
#  - Requer bastante RAM (a matriz densa intermediária tem ~55.737 genes x
#    16.291 células ≈ 7-8 GB antes de converter para esparsa;
#    recomenda-se >=16 GB de RAM livre).
# =============================================================================

library(data.table)
library(Matrix)
library(dplyr)
library(readr)
library(ggplot2)

set.seed(42)

raw_dir  <- here::here("01-melanoma-single-cell-biomarkers", "data", "raw")
proc_dir <- here::here("01-melanoma-single-cell-biomarkers", "data", "processed")
fig_dir  <- here::here("01-melanoma-single-cell-biomarkers", "figures")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

tpm_path <- file.path(raw_dir, "GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt")

# -----------------------------------------------------------------------------
# Carrega a matriz de expressão (data.table::fread é bem mais rápido/leve
# que read.csv/read_tsv para um arquivo deste tamanho)
# -----------------------------------------------------------------------------
message("Lendo matriz TPM (~4.5 GB) — isso pode levar vários minutos...")

header_line <- readLines(tpm_path, n = 1)
cell_ids <- strsplit(header_line, "\t")[[1]][-1]

expr_dt <- fread(tpm_path, skip = 2, header = FALSE, sep = "\t")
gene_names <- expr_dt[[1]]
expr_mat <- as.matrix(expr_dt[, -1, with = FALSE])
rownames(expr_mat) <- gene_names
colnames(expr_mat) <- cell_ids
rm(expr_dt)
gc()

expr_mat <- Matrix(expr_mat, sparse = TRUE)
message(sprintf("Matriz carregada: %d genes x %d células", nrow(expr_mat), ncol(expr_mat)))

# -----------------------------------------------------------------------------
# Reverte a transformação log2[(TPM/10)+1] -> TPM linear, só nos valores
# não-zero armazenados (mantém a matriz esparsa; evita materializar uma
# matriz densa de 2^x)
# -----------------------------------------------------------------------------
expr_tpm <- expr_mat
expr_tpm@x <- (2^(expr_tpm@x) - 1) * 10

# -----------------------------------------------------------------------------
# Métricas de QC por célula
# -----------------------------------------------------------------------------
qc_metrics <- tibble(
  cell_id      = colnames(expr_mat),
  n_genes      = Matrix::colSums(expr_mat > 0),
  total_tpm    = Matrix::colSums(expr_tpm)
)

mito_genes <- grep("^MT-", gene_names, value = TRUE)
if (length(mito_genes) > 0) {
  qc_metrics$pct_mito <- Matrix::colSums(expr_tpm[mito_genes, , drop = FALSE]) /
    qc_metrics$total_tpm * 100
} else {
  qc_metrics$pct_mito <- NA_real_
}

# -----------------------------------------------------------------------------
# Junta com os metadados de célula (response, patient_id, therapy)
# -----------------------------------------------------------------------------
cell_meta <- read_csv(file.path(proc_dir, "cell_metadata.csv"))
qc_metrics <- qc_metrics %>% left_join(cell_meta, by = "cell_id")

n_unmatched <- sum(is.na(qc_metrics$response))
if (n_unmatched > 0) {
  message(sprintf(
    "Aviso: %d células da matriz não encontraram metadado correspondente (response = NA).",
    n_unmatched
  ))
}

# -----------------------------------------------------------------------------
# Filtros de QC (thresholds iniciais — revisar após inspeção visual)
# -----------------------------------------------------------------------------
MIN_GENES <- 500
MAX_PCT_MITO <- 20

qc_metrics <- qc_metrics %>%
  mutate(
    pass_qc = n_genes >= MIN_GENES &
      (is.na(pct_mito) | pct_mito <= MAX_PCT_MITO) &
      !is.na(response)
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
  labs(title = "Genes detectados por célula (TPM > 0)", x = "n_genes", y = "count") +
  theme_minimal()
ggsave(file.path(fig_dir, "qc_n_genes.png"), p_genes, width = 6, height = 4)

p_mito <- ggplot(qc_metrics, aes(x = pct_mito, fill = pass_qc)) +
  geom_histogram(bins = 60) +
  geom_vline(xintercept = MAX_PCT_MITO, linetype = "dashed") +
  labs(title = "% mitocondrial por célula (TPM revertido)", x = "% mito", y = "count") +
  theme_minimal()
ggsave(file.path(fig_dir, "qc_pct_mito.png"), p_mito, width = 6, height = 4)

# -----------------------------------------------------------------------------
# Salva matriz filtrada (mantém a matriz log-transformada original, que é
# a forma esperada para o pré-processamento Seurat) + metadados
# -----------------------------------------------------------------------------
cells_keep <- qc_metrics %>% filter(pass_qc) %>% pull(cell_id)
log_tpm_filtered <- expr_mat[, cells_keep]

saveRDS(log_tpm_filtered, file.path(proc_dir, "log_tpm_filtered.rds"))
write_csv(qc_metrics, file.path(proc_dir, "qc_metrics.csv"))

message("QC concluído. Matriz filtrada salva em: ", proc_dir)
