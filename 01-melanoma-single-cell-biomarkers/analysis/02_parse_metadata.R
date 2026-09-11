# =============================================================================
# 02_parse_metadata.R
#
# GSE120575_patient_ID_single_cells.txt não é um CSV simples — é um
# arquivo no formato de template de submissão GEO (MINiML), com um bloco
# "SERIES" no topo, um bloco "SAMPLES" com uma linha por célula (a partir
# da linha que começa com "Sample name"), e blocos de protocolo/arquivo
# depois. Este script extrai só o bloco de dados por célula.
#
# Colunas relevantes no bloco de dados (por posição, já que o cabeçalho
# tem nomes longos/com erro de digitação no arquivo original):
#   X2 = title (= ID da célula, bate com o header da matriz TPM)
#   X5 = characteristics: patient ID (Pre=baseline; Post=on treatment)
#   X6 = characteristics: response (Responder / Non-responder)
#   X7 = characteristics: therapy
# =============================================================================

library(readr)
library(dplyr)
library(stringr)

set.seed(42)

raw_dir  <- here::here("01-melanoma-single-cell-biomarkers", "data", "raw")
proc_dir <- here::here("01-melanoma-single-cell-biomarkers", "data", "processed")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)

path <- file.path(raw_dir, "GSE120575_patient_ID_single_cells.txt")
lines <- read_lines(path)

# -----------------------------------------------------------------------------
# Localiza o bloco de dados por célula dinamicamente
# -----------------------------------------------------------------------------
header_idx <- which(str_starts(lines, "Sample name"))[1]
stopifnot(!is.na(header_idx))

data_start <- header_idx + 1
blank_rel <- which(lines[data_start:length(lines)] == "")[1]
data_end <- if (is.na(blank_rel)) length(lines) else (data_start + blank_rel - 2)

message(sprintf("Bloco de dados por célula: linhas %d a %d (%d células)",
                 data_start, data_end, data_end - data_start + 1))

data_block <- paste(lines[data_start:data_end], collapse = "\n")

meta_raw <- read_tsv(
  I(data_block),
  col_names = FALSE,
  col_types = cols(.default = "c")
)

# -----------------------------------------------------------------------------
# Seleciona e nomeia as colunas relevantes
# -----------------------------------------------------------------------------
cell_meta <- meta_raw %>%
  transmute(
    cell_id    = X2,
    patient_id = X5,
    response   = X6,
    therapy    = X7,
    timepoint  = case_when(
      str_starts(patient_id, "Pre_")  ~ "pre",
      str_starts(patient_id, "Post_") ~ "post",
      TRUE ~ NA_character_
    )
  )

stopifnot(!anyNA(cell_meta$cell_id))
message(sprintf("Metadados extraídos para %d células.", nrow(cell_meta)))
message("Distribuição de resposta:")
print(table(cell_meta$response, useNA = "ifany"))

write_csv(cell_meta, file.path(proc_dir, "cell_metadata.csv"))
message("Metadados salvos em: ", file.path(proc_dir, "cell_metadata.csv"))
