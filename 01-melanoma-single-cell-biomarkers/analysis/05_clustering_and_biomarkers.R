# =============================================================================
# 05_clustering_and_biomarkers.R
#
# Clustering (Louvain via Seurat), UMAP, e expressão diferencial entre
# Responders (R) e Non-responders (NR) — globalmente e por cluster.
#
# Nota: este dataset não vem com anotação de tipo celular publicada (ao
# contrário do GSE115978 usado no Projeto 02) — os clusters não
# supervisionados definidos aqui são a própria unidade de análise,
# seguindo a abordagem do paper original (definição de estados de
# células T por clustering).
# =============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(readr)
library(tibble)

set.seed(42)

proc_dir <- here::here("01-melanoma-single-cell-biomarkers", "data", "processed")
fig_dir  <- here::here("01-melanoma-single-cell-biomarkers", "figures")

seu <- readRDS(file.path(proc_dir, "seurat_preprocessed.rds"))

# -----------------------------------------------------------------------------
# Escolha do número de PCs (inspecionar ElbowPlot antes de fixar)
# -----------------------------------------------------------------------------
N_PCS <- 30

p_elbow <- ElbowPlot(seu, ndims = 50)
ggsave(file.path(fig_dir, "elbow_plot.png"), p_elbow, width = 6, height = 4)

# -----------------------------------------------------------------------------
# Clustering + UMAP
# -----------------------------------------------------------------------------
seu <- FindNeighbors(seu, dims = 1:N_PCS)
seu <- FindClusters(seu, resolution = 0.8, random.seed = 42)
seu <- RunUMAP(seu, dims = 1:N_PCS, seed.use = 42)

p_umap_cluster <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters",
                           label = TRUE) +
  labs(title = "UMAP — clusters não supervisionados")
ggsave(file.path(fig_dir, "umap_clusters.png"), p_umap_cluster, width = 7, height = 6)

p_umap_response <- DimPlot(seu, reduction = "umap", group.by = "response") +
  labs(title = "UMAP — resposta clínica (R/NR)")
ggsave(file.path(fig_dir, "umap_response.png"), p_umap_response, width = 7, height = 6)

if ("therapy" %in% colnames(seu@meta.data)) {
  p_umap_therapy <- DimPlot(seu, reduction = "umap", group.by = "therapy") +
    labs(title = "UMAP — linha de terapia")
  ggsave(file.path(fig_dir, "umap_therapy.png"), p_umap_therapy, width = 7, height = 6)
}

# -----------------------------------------------------------------------------
# Composição de cada cluster por grupo de resposta
# -----------------------------------------------------------------------------
crosstab <- table(cluster = seu@meta.data$seurat_clusters, response = seu@meta.data$response)
write.csv(as.data.frame.matrix(crosstab), file.path(proc_dir, "cluster_vs_response_crosstab.csv"))

saveRDS(seu, file.path(proc_dir, "seurat_clustered.rds"))

# -----------------------------------------------------------------------------
# DE global R vs NR (todas as células)
# -----------------------------------------------------------------------------
Idents(seu) <- seu@meta.data$response

de_global <- FindMarkers(
  seu,
  ident.1 = "Responder",
  ident.2 = "Non-responder",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1
) %>%
  rownames_to_column("gene") %>%
  arrange(p_val_adj)

write_csv(de_global, file.path(proc_dir, "de_global_R_vs_NR.csv"))

# -----------------------------------------------------------------------------
# DE por cluster (R vs NR dentro de cada cluster não supervisionado)
# -----------------------------------------------------------------------------
clusters <- unique(seu@meta.data$seurat_clusters)
de_by_cluster <- list()

for (cl in clusters) {
  sub <- subset(seu, subset = seurat_clusters == cl)
  Idents(sub) <- sub@meta.data$response

  groups_present <- table(Idents(sub))
  if (length(groups_present) < 2 || any(groups_present < 3)) next

  de_cl <- tryCatch(
    FindMarkers(
      sub,
      ident.1 = "Responder",
      ident.2 = "Non-responder",
      test.use = "wilcox",
      logfc.threshold = 0.25,
      min.pct = 0.1
    ) %>%
      rownames_to_column("gene") %>%
      mutate(cluster = cl) %>%
      arrange(p_val_adj),
    error = function(e) NULL
  )
  if (!is.null(de_cl)) de_by_cluster[[as.character(cl)]] <- de_cl
}

de_by_cluster_df <- bind_rows(de_by_cluster)
write_csv(de_by_cluster_df, file.path(proc_dir, "de_by_cluster_R_vs_NR.csv"))

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

message("Clustering e análise de biomarcadores concluídos.")
message("Top 10 genes candidatos (global R vs NR):")
print(head(de_global %>% filter(sig != "NS"), 10))
