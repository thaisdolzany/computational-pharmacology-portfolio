# =============================================================================
# 04_clustering_annotation.R
#
# Graph-based clustering (Louvain), UMAP, and comparison with the published
# cell-type labels (column cell.types in GSE115978_cell.annotations.csv).
# =============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)

set.seed(42)

proc_dir <- here::here("melanoma-scrnaseq-resistance", "data", "processed")
fig_dir  <- here::here("melanoma-scrnaseq-resistance", "figures")

seu <- readRDS(file.path(proc_dir, "seurat_preprocessed.rds"))

# 30 PCs, chosen from the elbow plot
N_PCS <- 30
ggsave(file.path(fig_dir, "elbow_plot.png"),
       ElbowPlot(seu, ndims = 50) + labs(title = "Variance explained by each PC"),
       width = 6, height = 4, bg = "white")

seu <- FindNeighbors(seu, dims = 1:N_PCS)
seu <- FindClusters(seu, resolution = 0.8, random.seed = 42)
seu <- RunUMAP(seu, dims = 1:N_PCS, seed.use = 42)

umap <- function(group, title, file) {
  p <- DimPlot(seu, reduction = "umap", group.by = group, label = group == "seurat_clusters") +
    labs(title = title, x = "UMAP 1", y = "UMAP 2")
  ggsave(file.path(fig_dir, file), p, width = 7, height = 6, bg = "white")
}
umap("seurat_clusters", "Unsupervised clusters", "umap_clusters.png")
umap("cell.types", "Published cell-type labels", "umap_celltypes.png")
umap("treatment.group", "Biopsy timing (treatment-naive or post-treatment)", "umap_treatment_group.png")
umap("Cohort", "Cohort of origin (Tirosh or New)", "umap_cohort.png")

crosstab <- table(cluster = seu$seurat_clusters, cell_type_published = seu$cell.types)
write.csv(as.data.frame.matrix(crosstab), file.path(proc_dir, "cluster_vs_celltype_crosstab.csv"))

saveRDS(seu, file.path(proc_dir, "seurat_clustered.rds"))
message(sprintf("Clustering done: %d clusters", length(levels(seu$seurat_clusters))))
