# =============================================================================
# 05_resistance_program_analysis.R
#
# Exploratory, cell-level analysis of malignant cells, following the logic of
# Jerby-Arnon et al. (Cell 2018) within what the public data allow:
#
#   (a) differential expression in malignant cells, post-treatment vs
#       treatment-naive biopsies (Wilcoxon, cells as observations)
#   (b) a per-sample signature score from the genes up in (a), correlated with
#       the fraction of T cells in the same sample
#
# Cells from the same tumour are not independent, so the p-values in (a) are
# inflated. 06_pseudobulk_de.R repeats the comparison with tumours as the unit
# of replication and is the analysis to use for inference.
#
# The dataset has no clinical response labels (responder / non-responder).
# =============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(readr)
library(tibble)

set.seed(42)

proc_dir <- here::here("melanoma-scrnaseq-resistance", "data", "processed")
fig_dir  <- here::here("melanoma-scrnaseq-resistance", "figures")

seu <- readRDS(file.path(proc_dir, "seurat_clustered.rds"))
stopifnot(all(c("cell.types", "treatment.group", "samples") %in% colnames(seu@meta.data)))

# -----------------------------------------------------------------------------
# (a) Malignant cells: post-treatment vs treatment-naive
# -----------------------------------------------------------------------------
mal <- subset(seu, subset = cell.types == "Mal")
Idents(mal) <- mal$treatment.group

de_mal <- FindMarkers(mal, ident.1 = "post.treatment", ident.2 = "treatment.naive",
                      test.use = "wilcox", logfc.threshold = 0.25, min.pct = 0.1) %>%
  rownames_to_column("gene") %>%
  arrange(p_val_adj)
write_csv(de_mal, file.path(proc_dir, "de_malignant_post_vs_naive.csv"))

p_volcano <- de_mal %>%
  mutate(group = case_when(p_val_adj < 0.05 & avg_log2FC > 0.25 ~ "Higher post-treatment",
                           p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Lower post-treatment",
                           TRUE ~ "Not significant"),
         neglog_padj = -log10(pmax(p_val_adj, 1e-300))) %>%
  ggplot(aes(x = avg_log2FC, y = neglog_padj, colour = group)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_colour_manual(values = c("Higher post-treatment" = "#D55E00",
                                 "Lower post-treatment" = "#0072B2",
                                 "Not significant" = "grey75")) +
  labs(title = "Malignant cells, post-treatment vs treatment-naive (cell-level)",
       subtitle = "Wilcoxon test with cells as observations; p-values inflated, see pseudobulk analysis",
       x = "Average log2 fold change", y = "-log10 adjusted p-value", colour = NULL) +
  theme_minimal()
ggsave(file.path(fig_dir, "volcano_malignant_post_vs_naive.png"), p_volcano,
       width = 7.5, height = 6, bg = "white")

# -----------------------------------------------------------------------------
# (b) Signature score per sample vs T-cell fraction in the same sample
# -----------------------------------------------------------------------------
resistance_genes <- de_mal %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.25) %>%
  slice_max(order_by = avg_log2FC, n = 50) %>%
  pull(gene)

if (length(resistance_genes) >= 5) {
  mal <- AddModuleScore(mal, features = list(resistance_genes), name = "resistance_score", seed = 42)

  score_by_sample <- mal@meta.data %>%
    group_by(samples) %>%
    summarise(mean_resistance_score = mean(resistance_score1), .groups = "drop")

  tcell_frac <- seu@meta.data %>%
    group_by(samples) %>%
    summarise(n_cells = n(),
              n_tcell = sum(grepl("^T\\.", cell.types)),
              tcell_fraction = n_tcell / n_cells,
              .groups = "drop")

  sample_summary <- inner_join(score_by_sample, tcell_frac, by = "samples")
  write_csv(sample_summary, file.path(proc_dir, "resistance_score_vs_tcell_fraction.csv"))

  ct <- cor.test(sample_summary$mean_resistance_score, sample_summary$tcell_fraction,
                 method = "spearman", exact = FALSE)
  message(sprintf("Spearman, score vs T-cell fraction: rho = %.3f, p = %.3f, n = %d samples",
                  ct$estimate, ct$p.value, nrow(sample_summary)))

  p_cor <- ggplot(sample_summary, aes(x = mean_resistance_score, y = tcell_fraction)) +
    geom_point(size = 2, colour = "#0072B2") +
    labs(title = "Malignant-cell signature score vs T-cell fraction, per sample",
         subtitle = sprintf("Spearman rho = %.2f, p = %.2f, n = %d samples",
                            ct$estimate, ct$p.value, nrow(sample_summary)),
         x = "Mean signature score in malignant cells", y = "Fraction of T cells in the sample") +
    theme_minimal()
  ggsave(file.path(fig_dir, "resistance_score_vs_tcell_fraction.png"), p_cor,
         width = 7, height = 5.5, bg = "white")
} else {
  message("Fewer than 5 significant genes: signature score not computed.")
}

message("Top 10 genes (cell-level, post-treatment vs treatment-naive):")
print(head(de_mal %>% filter(p_val_adj < 0.05), 10))
