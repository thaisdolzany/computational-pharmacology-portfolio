# =============================================================================
# 05_biomarker_analysis.R
#
# Expressão diferencial entre responders (R) e non-responders (NR) por
# população celular, buscando candidatos a biomarcador preditivo de
# resposta a checkpoint inhibitors.
# =============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(readr)

set.seed(42)

proc_dir <- here::here("01-melanoma-single-cell-biomarkers", "data", "processed")
fig_dir  <- here::here("01-melanoma-single-cell-biomarkers", "figures")

seu <- readRDS(file.path(proc_dir, "seurat_clustered.rds"))

stopifnot("treatment.group" %in% colnames(seu@meta.data))
Idents(seu) <- seu@meta.data$treatment.group

# -----------------------------------------------------------------------------
# DE global R vs NR (todas as células) — teste Wilcoxon (default do Seurat)
# -----------------------------------------------------------------------------
de_global <- FindMarkers(
  seu,
  ident.1 = "Responder",
  ident.2 = "Non-responder",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1
) %>%
  tibble::rownames_to_column("gene") %>%
  arrange(p_val_adj)

write_csv(de_global, file.path(proc_dir, "de_global_R_vs_NR.csv"))

# -----------------------------------------------------------------------------
# DE por população celular (subtipo publicado), quando disponível
# -----------------------------------------------------------------------------
if ("cell.types" %in% colnames(seu@meta.data)) {
  cell_types <- unique(seu@meta.data$cell.types)
  de_by_celltype <- list()

  for (ct in cell_types) {
    sub <- subset(seu, subset = cell.types == ct)
    Idents(sub) <- sub@meta.data$treatment.group

    groups_present <- table(Idents(sub))
    if (length(groups_present) < 2 || any(groups_present < 3)) next

    de_ct <- tryCatch(
      FindMarkers(
        sub,
        ident.1 = "Responder",
        ident.2 = "Non-responder",
        test.use = "wilcox",
        logfc.threshold = 0.25,
        min.pct = 0.1
      ) %>%
        tibble::rownames_to_column("gene") %>%
        mutate(cell_type = ct) %>%
        arrange(p_val_adj),
      error = function(e) NULL
    )
    if (!is.null(de_ct)) de_by_celltype[[ct]] <- de_ct
  }

  de_by_celltype_df <- bind_rows(de_by_celltype)
  write_csv(de_by_celltype_df,
            file.path(proc_dir, "de_by_celltype_R_vs_NR.csv"))
}

# -----------------------------------------------------------------------------
# Volcano plot global
# -----------------------------------------------------------------------------
de_global <- de_global %>%
  mutate(
    sig = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.25  ~ "Up em Responders",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Up em Non-responders",
      TRUE ~ "NS"
    )
  )

p_volcano <- ggplot(de_global, aes(x = avg_log2FC, y = -log10(p_val_adj), color = sig)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c(
    "Up em Responders" = "#2166AC",
    "Up em Non-responders" = "#B2182B",
    "NS" = "grey70"
  )) +
  labs(
    title = "Expressão diferencial: Responders vs Non-responders",
    x = "avg log2FC", y = "-log10(p_val_adj)", color = NULL
  ) +
  theme_minimal()
ggsave(file.path(fig_dir, "volcano_R_vs_NR.png"), p_volcano, width = 7, height = 6)

message("Análise de biomarcadores concluída.")
message("Top 10 genes candidatos:")
print(head(de_global %>% filter(sig != "NS"), 10))
