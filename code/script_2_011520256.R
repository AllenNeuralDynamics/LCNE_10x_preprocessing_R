
rm(list = ls())
gc()  # Garbage collection at start

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)


plot_integration_results <- function(combined) {
  p1 <- DimPlot(combined, reduction = "pca", group.by = "orig.ident", pt.size = 0.1) +
    ggtitle("Before Harmony Integration")
  p2 <- DimPlot(combined, reduction = "harmony", group.by = "orig.ident", pt.size = 0.1) +
    ggtitle("After Harmony Integration")
  p3 <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 0.1) +
    ggtitle("Clusters")
  return(list(batch_effects = p1 + p2, clusters = p3))
}



# load the results 
combined <- readRDS("../scratch/LC_clustered_harmony.rds")

## not neeeded. 
# integration_plots <- plot_integration_results(combined)
# print(integration_plots$batch_effects)

print(DimPlot(combined, reduction = "umap", group.by = "seurat_clusters",
              label = TRUE, repel = TRUE, pt.size = 0.3))
DimPlot(combined, reduction = "umap", group.by = "batch", pt.size = 0.3)




lc_markers <- c("Dbh", "Slc6a2", "Th", "Slc18a2", "Tacr3")
FeaturePlot( combined, features = lc_markers,
  reduction = "umap", pt.size = 0.3, ncol = 3)

avg_expr <- AverageExpression(combined,features = lc_markers,
  assays = "RNA",  slot = "data")$RNA



DotPlot(
  combined,
  features = lc_markers,
  group.by = "seurat_clusters"
) + RotatedAxis()



clusters_of_interest <- c("9", "26")
VlnPlot(combined,
  features = lc_markers,
  idents = clusters_of_interest,
  pt.size = 0,ncol = 3)


Idents(combined) <- "seurat_clusters"
cluster_interest <- c(9)  # or c(9, 24).   # 5151 samples slected
combined_lc <- subset(combined, idents = cluster_interest)

dim(combined_lc) # 5151





saveRDS(combined_lc, "../scratch/LC_cluster_interest.rds")




