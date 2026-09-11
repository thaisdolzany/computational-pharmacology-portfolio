# =============================================================================
# 01_download_gse120575.R
#
# Baixa os arquivos supplementary do GEO GSE120575 (Sade-Feldman et al.,
# Cell 2018, PMID 30388456) para data/raw/.
#
# Dataset: scRNA-seq (Smart-seq2) de 16.291 células CD45+ de 48 biópsias
# de melanoma (32 pacientes) tratados com checkpoint inhibitors, com
# rótulo de resposta clínica (Responder/Non-responder) por célula.
#
# AVISO: a matriz TPM descompactada tem ~4,5 GB. Garanta espaço em disco
# livre (>=10 GB recomendado, considerando o .gz + o .txt descompactado).
# =============================================================================

library(GEOquery)

set.seed(42)

raw_dir <- here::here("01-melanoma-single-cell-biomarkers", "data", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

geo_acc <- "GSE120575"

# GEOquery baixa os arquivos supplementary associados ao accession:
#   GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz
#   GSE120575_patient_ID_single_cells.txt.gz
getGEOSuppFiles(
  GEO = geo_acc,
  baseDir = raw_dir,
  makeDirectory = FALSE,
  fetch_files = TRUE
)

# Descompacta os .gz para leitura direta nos próximos scripts
gz_files <- list.files(raw_dir, pattern = "\\.gz$", full.names = TRUE)
lapply(gz_files, R.utils::gunzip, remove = FALSE, overwrite = TRUE)

message("Download concluído. Arquivos em: ", raw_dir)
message("Arquivos esperados:")
message("  - ", geo_acc, "_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt (~4.5 GB)")
message("  - ", geo_acc, "_patient_ID_single_cells.txt")
