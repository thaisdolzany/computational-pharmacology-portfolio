# =============================================================================
# 01_download_gse115978.R
#
# Baixa os arquivos supplementary do GEO GSE115978 (Jerby-Arnon et al.,
# Cell 2018, PMID 30388455) para data/raw/.
#
# Dataset: scRNA-seq (Smart-seq2) de 31 tumores de melanoma, combinando as
# coortes Tirosh (GSE72056) e New. Usado para identificar um programa de
# resistência em células malignas associado a exclusão de células T.
# =============================================================================

library(GEOquery)

set.seed(42)

raw_dir <- here::here("02-melanoma-drug-resistance", "data", "raw")
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
