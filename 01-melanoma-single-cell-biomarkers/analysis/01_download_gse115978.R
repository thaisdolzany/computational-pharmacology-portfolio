# =============================================================================
# 01_download_gse115978.R
#
# Baixa os arquivos supplementary do GEO GSE115978 (Sade-Feldman et al.,
# Cell 2018) para data/raw/.
#
# Dataset: scRNA-seq (Smart-seq2) de melanoma metastático, pré/pós
# checkpoint inhibitor, com anotação de resposta clínica (R/NR) por célula.
# =============================================================================

library(GEOquery)

set.seed(42)

raw_dir <- here::here("01-melanoma-single-cell-biomarkers", "data", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

geo_acc <- "GSE115978"

# GEOquery baixa todos os arquivos supplementary associados ao accession.
# Inclui: GSE115978_counts.csv.gz, GSE115978_tpm.csv.gz,
#         GSE115978_cell.annotations.csv.gz
getGEOSuppFiles(
  GEO = geo_acc,
  baseDir = raw_dir,
  makeDirectory = FALSE,
  fetch_files = TRUE
)

# Descompacta os .csv.gz para leitura direta nos próximos scripts
gz_files <- list.files(raw_dir, pattern = "\\.gz$", full.names = TRUE)
lapply(gz_files, R.utils::gunzip, remove = FALSE, overwrite = TRUE)

message("Download concluído. Arquivos em: ", raw_dir)
message("Arquivos esperados:")
message("  - ", geo_acc, "_counts.csv")
message("  - ", geo_acc, "_tpm.csv")
message("  - ", geo_acc, "_cell.annotations.csv")
