# =============================================================================
# 01_download_gse115978.R
#
# Downloads the GEO supplementary files of GSE115978 (Jerby-Arnon et al.,
# Cell 2018, PMID 30388455) into data/raw/.
#
# Dataset: Smart-seq2 scRNA-seq of 31 melanoma tumours, combining the Tirosh
# cohort (GSE72056) with new samples from the same study.
# =============================================================================

library(GEOquery)

raw_dir <- here::here("melanoma-scrnaseq-resistance", "data", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

geo_acc <- "GSE115978"

# Supplementary files: GSE115978_counts.csv.gz, GSE115978_tpm.csv.gz,
# GSE115978_cell.annotations.csv.gz
getGEOSuppFiles(
  GEO = geo_acc,
  baseDir = raw_dir,
  makeDirectory = FALSE,
  fetch_files = TRUE
)

gz_files <- list.files(raw_dir, pattern = "\\.gz$", full.names = TRUE)
lapply(gz_files, R.utils::gunzip, remove = FALSE, overwrite = TRUE)

message("Download complete: ", raw_dir)
