# =============================================================================
# 04_clustering_annotation.R
#
# Clustering (Louvain via Seurat), UMAP, e validação cruzada com as
# anotações de tipo celular publicadas pelos autores (coluna cell.types
# em GSE115978_cell.annotations.csv).
# =============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)

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

# -----------------------------------------------------------------------------
# Visualizações
# -----------------------------------------------------------------------------
p_umap_cluster <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters",
                           label = TRUE) +
  labs(title = "UMAP — clusters não supervisionados")
ggsave(file.path(fig_dir, "umap_clusters.png"), p_umap_cluster, width = 7, height = 6)

if ("cell.types" %in% colnames(seu@meta.data)) {
  p_umap_celltype <- DimPlot(seu, reduction = "umap", group.by = "cell.types") +
    labs(title = "UMAP — anotação de tipo celular (publicada)")
  ggsave(file.path(fig_dir, "umap_celltypes.png"), p_umap_celltype, width = 7, height = 6)
}

if ("treatment.group" %in% colnames(seu@meta.data)) {
  p_umap_response <- DimPlot(seu, reduction = "umap", group.by = "treatment.group") +
    labs(title = "UMAP — grupo de resposta clínica (R/NR)")
  ggsave(file.path(fig_dir, "umap_response.png"), p_umap_response, width = 7, height = 6)
}

# -----------------------------------------------------------------------------
# Tabela de correspondência cluster x tipo celular publicado
# -----------------------------------------------------------------------------
if ("cell.types" %in% colnames(seu@meta.data)) {
  crosstab <- table(
    cluster = seu@meta.data$seurat_clusters,
    cell_type_published = seu@meta.data$cell.types
  )
  write.csv(as.data.frame.matrix(crosstab),
            file.path(proc_dir, "cluster_vs_celltype_crosstab.csv"))
}

saveRDS(seu, file.path(proc_dir, "seurat_clustered.rds"))

message("Clustering e anotação concluídos. Objeto salvo em: ", proc_dir)
