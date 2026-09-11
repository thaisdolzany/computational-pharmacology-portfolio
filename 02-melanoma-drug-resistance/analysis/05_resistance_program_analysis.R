# =============================================================================
# 05_resistance_program_analysis.R
#
# Investiga um programa transcricional em células malignas associado a
# exclusão de células T, seguindo a lógica do paper original (Jerby-Arnon
# et al., Cell 2018), mas simplificado ao que os dados públicos permitem:
#
#   (a) DE em células malignas: post.treatment vs treatment.naive
#       (proxy de programa selecionado/induzido por tratamento)
#   (b) Score de assinatura por amostra (malignant cells) correlacionado
#       com a fração de células T na mesma amostra (proxy de exclusão
#       imune — quanto maior o score, menor a infiltração de T esperada
#       se o programa de fato promove exclusão)
#
# Nota: este dataset não tem rótulo de resposta clínica (R/NR) — ver
# README, seção "Dataset" e "Limitations".
# =============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(readr)
library(tibble)

set.seed(42)

proc_dir <- here::here("02-melanoma-drug-resistance", "data", "processed")
fig_dir  <- here::here("02-melanoma-drug-resistance", "figures")

seu <- readRDS(file.path(proc_dir, "seurat_clustered.rds"))

stopifnot(all(c("cell.types", "treatment.group", "samples") %in% colnames(seu@meta.data)))

# -----------------------------------------------------------------------------
# (a) DE em células malignas: post.treatment vs treatment.naive
# -----------------------------------------------------------------------------
mal <- subset(seu, subset = cell.types == "Mal")
Idents(mal) <- mal@meta.data$treatment.group

de_mal <- FindMarkers(
  mal,
  ident.1 = "post.treatment",
  ident.2 = "treatment.naive",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1
) %>%
  rownames_to_column("gene") %>%
  arrange(p_val_adj)

write_csv(de_mal, file.path(proc_dir, "de_malignant_post_vs_naive.csv"))

p_volcano <- de_mal %>%
  mutate(sig = p_val_adj < 0.05 & abs(avg_log2FC) > 0.25) %>%
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), color = sig)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("TRUE" = "#B2182B", "FALSE" = "grey70")) +
  labs(
    title = "Células malignas: post.treatment vs treatment.naive",
    x = "avg log2FC", y = "-log10(p_val_adj)", color = "Significativo"
  ) +
  theme_minimal()
ggsave(file.path(fig_dir, "volcano_malignant_post_vs_naive.png"), p_volcano,
       width = 7, height = 6)

# -----------------------------------------------------------------------------
# (b) Score de assinatura de resistência por célula maligna, agregado por
#     amostra, correlacionado com fração de células T na amostra
# -----------------------------------------------------------------------------
resistance_genes <- de_mal %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.25) %>%
  slice_max(order_by = avg_log2FC, n = 50) %>%
  pull(gene)

if (length(resistance_genes) >= 5) {
  mal <- AddModuleScore(
    mal,
    features = list(resistance_genes),
    name = "resistance_score",
    seed = 42
  )

  score_by_sample <- mal@meta.data %>%
    group_by(samples) %>%
    summarise(mean_resistance_score = mean(resistance_score1), .groups = "drop")

  # Fração de células T (T.CD8 + T.CD4, se presentes) por amostra, no objeto
  # completo (todas as células, não só malignas)
  tcell_frac <- seu@meta.data %>%
    mutate(is_tcell = grepl("^T\\.", cell.types)) %>%
    group_by(samples) %>%
    summarise(
      n_cells = n(),
      n_tcell = sum(is_tcell),
      tcell_fraction = n_tcell / n_cells,
      .groups = "drop"
    )

  sample_summary <- inner_join(score_by_sample, tcell_frac, by = "samples")
  write_csv(sample_summary, file.path(proc_dir, "resistance_score_vs_tcell_fraction.csv"))

  cor_test <- cor.test(
    sample_summary$mean_resistance_score,
    sample_summary$tcell_fraction,
    method = "spearman"
  )
  message(sprintf(
    "Correlação (Spearman) score de resistência x fração de células T: rho = %.3f, p = %.4f",
    cor_test$estimate, cor_test$p.value
  ))

  p_cor <- ggplot(sample_summary, aes(x = mean_resistance_score, y = tcell_fraction)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "#2166AC") +
    labs(
      title = "Score de programa de resistência vs fração de células T por amostra",
      x = "Score médio de resistência (células malignas)",
      y = "Fração de células T na amostra"
    ) +
    theme_minimal()
  ggsave(file.path(fig_dir, "resistance_score_vs_tcell_fraction.png"), p_cor,
         width = 7, height = 6)
} else {
  message("Menos de 5 genes significativos no DE — score de assinatura não calculado.")
}

message("Análise de programa de resistência concluída.")
message("Top 10 genes candidatos (post.treatment vs treatment.naive):")
print(head(de_mal %>% filter(p_val_adj < 0.05), 10))
