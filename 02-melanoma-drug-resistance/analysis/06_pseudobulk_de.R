# =============================================================================
# 06_pseudobulk_de.R
#
# Refaz o DE de células malignas (post.treatment vs treatment.naive) no nível
# da amostra, que é a unidade real de replicação.
#
# Por que: o script 05 usa Wilcoxon célula a célula (FindMarkers). O tratamento
# é uma característica do tumor, não da célula, e células malignas de melanoma
# agrupam fortemente por paciente. Tratar milhares de células como réplicas
# independentes infla os p-valores (pseudo-replicação) e deixa genes
# específicos de um ou dois pacientes parecerem efeitos de tratamento.
#
# Abordagem: soma das contagens brutas das células malignas por amostra
# (pseudobulk), DESeq2 com desenho ~ Cohort + treatment.group (Cohort controla
# o lote Tirosh vs New, que não está confundido com o tratamento).
# Só entram amostras com >= MIN_CELLS células malignas.
# =============================================================================

library(Matrix)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(readr)
library(tibble)
library(tidyr)

set.seed(42)

MIN_CELLS <- 10
GENES_05  <- c("MIA", "EMP1", "SERPINF1", "PFN1")  # destaques do README (script 05)

proc_dir <- here::here("02-melanoma-drug-resistance", "data", "processed")
raw_dir  <- here::here("02-melanoma-drug-resistance", "data", "raw")
res_dir  <- here::here("02-melanoma-drug-resistance", "results")
fig_dir  <- here::here("02-melanoma-drug-resistance", "figures")
dir.create(res_dir, showWarnings = FALSE)

counts <- readRDS(file.path(proc_dir, "counts_filtered.rds"))
ann    <- read_csv(file.path(raw_dir, "GSE115978_cell.annotations.csv"), show_col_types = FALSE)

mal <- ann %>% filter(cell.types == "Mal", cells %in% colnames(counts))
stopifnot(all(mal$cells %in% colnames(counts)))
stopifnot(all(as.matrix(counts[1:100, 1:100]) == round(as.matrix(counts[1:100, 1:100]))))

# -----------------------------------------------------------------------------
# Pseudobulk: soma das contagens por amostra
# -----------------------------------------------------------------------------
sample_info <- mal %>%
  count(samples, treatment.group, Cohort, name = "n_malignant_cells") %>%
  filter(n_malignant_cells >= MIN_CELLS) %>%
  mutate(
    treatment.group = factor(treatment.group, levels = c("treatment.naive", "post.treatment")),
    Cohort = factor(Cohort)
  )
mal <- mal %>% filter(samples %in% sample_info$samples)

membership <- sparseMatrix(
  i = match(mal$cells, colnames(counts)),
  j = match(mal$samples, sample_info$samples),
  x = 1,
  dims = c(ncol(counts), nrow(sample_info))
)
pb <- as.matrix(counts %*% membership)
colnames(pb) <- sample_info$samples
storage.mode(pb) <- "integer"

message(sprintf("Pseudobulk: %d genes x %d amostras (%s)", nrow(pb), ncol(pb),
                paste(names(table(sample_info$treatment.group)), table(sample_info$treatment.group),
                      sep = " = ", collapse = ", ")))

# filtro de genes pouco expressos: >= 10 contagens em pelo menos n_min amostras,
# onde n_min = tamanho do menor grupo
n_min <- min(table(sample_info$treatment.group))
keep  <- rowSums(pb >= 10) >= n_min
pb    <- pb[keep, ]
message(sprintf("Genes após filtro: %d", nrow(pb)))

# -----------------------------------------------------------------------------
# DESeq2
# -----------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = pb,
  colData   = column_to_rownames(as.data.frame(sample_info), "samples")[colnames(pb), ],
  design    = ~ Cohort + treatment.group
)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("treatment.group", "post.treatment", "treatment.naive"), alpha = 0.05)

de_pb <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  arrange(padj)
write_csv(de_pb, file.path(res_dir, "pseudobulk_de_malignant_post_vs_naive.csv"))
write_csv(sample_info, file.path(res_dir, "pseudobulk_samples.csv"))

# -----------------------------------------------------------------------------
# Comparação com o DE célula a célula (script 05)
# -----------------------------------------------------------------------------
de_cell <- read_csv(file.path(proc_dir, "de_malignant_post_vs_naive.csv"), show_col_types = FALSE)
sig_cell <- de_cell %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>% pull(gene)
sig_pb   <- de_pb %>% filter(!is.na(padj), padj < 0.05) %>% pull(gene)

comparison <- tibble(
  n_sig_cell_level        = length(sig_cell),
  n_sig_pseudobulk        = length(sig_pb),
  n_sig_both              = length(intersect(sig_cell, sig_pb)),
  pct_cell_hits_confirmed = round(100 * length(intersect(sig_cell, sig_pb)) / length(sig_cell), 1)
)
write_csv(comparison, file.path(res_dir, "cell_vs_pseudobulk_comparison.csv"))
print(comparison)

highlight <- de_pb %>%
  filter(gene %in% GENES_05) %>%
  select(gene, baseMean, log2FoldChange, lfcSE, pvalue, padj) %>%
  left_join(de_cell %>% select(gene, avg_log2FC, p_val_adj), by = "gene")
write_csv(highlight, file.path(res_dir, "pseudobulk_vs_cell_highlight_genes.csv"))
print(highlight)

message("Top 15 genes pseudobulk:")
print(head(de_pb %>% select(gene, log2FoldChange, padj), 15))

# -----------------------------------------------------------------------------
# Figura: expressão por tumor dos genes destacados no script 05
# (cada ponto = uma amostra; é a unidade de replicação)
# -----------------------------------------------------------------------------
vsd <- vst(dds, blind = FALSE)
genes_plot <- intersect(GENES_05, rownames(vsd))
expr <- as.data.frame(t(assay(vsd)[genes_plot, , drop = FALSE])) %>%
  rownames_to_column("samples") %>%
  pivot_longer(-samples, names_to = "gene", values_to = "vst") %>%
  left_join(sample_info, by = "samples") %>%
  left_join(highlight %>% select(gene, padj), by = "gene") %>%
  mutate(label = sprintf("%s\npseudobulk padj = %s", gene, signif(padj, 2)),
         group = recode(treatment.group, treatment.naive = "Naive", post.treatment = "Post-treatment"))

p <- ggplot(expr, aes(x = group, y = vst)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, colour = "grey40") +
  geom_point(aes(shape = Cohort, colour = group), position = position_jitter(width = 0.12, seed = 42), size = 2) +
  scale_colour_manual(values = c("Naive" = "#0072B2", "Post-treatment" = "#D55E00"), guide = "none") +
  facet_wrap(~ label, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = "Expression (DESeq2 VST, pseudobulk per tumour)",
       title = "Malignant cells, one point per tumour",
       subtitle = sprintf("n = %d naive, %d post-treatment tumours (>= %d malignant cells each)",
                          sum(sample_info$treatment.group == "treatment.naive"),
                          sum(sample_info$treatment.group == "post.treatment"), MIN_CELLS)) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.x = element_blank())
ggsave(file.path(fig_dir, "pseudobulk_highlight_genes.png"), p, width = 9, height = 3.8, dpi = 300, bg = "white")

writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_06.txt"))
message("Pseudobulk DE concluído.")
