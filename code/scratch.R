# jsut the scratch to combine 03 subclustering with their corresponding plotting functions. 
rm(list = ls())

input_file = "../scratch/LC_cluster_interest.rds"
output_file <- "../scratch/LC_subclusters_filtered.rds"

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
source("./func/integration_functions.R")

combined_lc <- readRDS(input_file)


# LC only preprocessing

combined_lc <- NormalizeData(combined_lc)
combined_lc <- FindVariableFeatures(combined_lc, nfeatures = 2000)
combined_lc <- ScaleData(combined_lc)

combined_lc <- RunPCA(combined_lc, npcs = 30)
combined_lc <- RunUMAP(combined_lc, dims = 1:10, return.model = TRUE)

combined_lc <- FindNeighbors(combined_lc, dims = 1:10)
combined_lc <- FindClusters(combined_lc, resolution = 0.4)



# QC filtering steps

combined_lc$nCount_RNA   <- as.numeric(combined_lc$nCount_RNA)
combined_lc$nFeature_RNA <- as.numeric(combined_lc$nFeature_RNA)
combined_lc$percent.mt   <- as.numeric(combined_lc$percent.mt)
combined_lc$doublet_score   <- as.numeric(combined_lc$doublet_score)

qc_df <- SeuratObject::FetchData(
  combined_lc,
  vars = c("nCount_RNA", "nFeature_RNA", "percent.mt", "seurat_clusters"),
  layer = "data"  # or whatever layer your QC metrics live in
)
qc_df$seurat_clusters <- as.factor(qc_df$seurat_clusters)


## plotting for combined_LC (the entire thing)
lc_markers <- c("Dbh", "Slc6a2", "Th", "Slc18a2", "Tacr3")

FeaturePlot(combined_lc,
  features  = lc_markers,
  reduction = "umap", pt.size   = 0.3,ncol      = 3)

# Also see the subclusters overlaid
DimPlot(combined_lc,
  reduction = "umap", group.by  = "seurat_clusters",
  label     = TRUE, repel = TRUE, pt.size = 0.4)






# some previous plotting with interested clusters -----------------------------------------------------


clusters_to_check <- c("0","2", "5", "7","8")

combined_lc$cluster_5_7_vs_other <- ifelse(
  combined_lc$seurat_clusters %in% c("5", "7"),
  "5_7_candidate",
  "other"
)




## QC for these clusters
VlnPlot(
  subset(combined_lc, idents = clusters_to_check),
  features = c("nCount_RNA", "nFeature_RNA", "doublet_score"),
  group.by = "seurat_clusters",
  pt.size  = 0.1,
  ncol     = 3
)




VlnPlot(
  subset(combined_lc, idents = clusters_to_check),
  features = lc_markers,
  group.by = "seurat_clusters",
  pt.size  = 0.1,
  ncol     = 3
)


DotPlot(
  combined_lc,
  features = lc_markers,
  group.by = "cluster_5_7_vs_other"
) + RotatedAxis()



bad_subclusters <- c("8","7")
combined_lc_filtered <- subset(combined_lc, idents = setdiff(levels(Idents(combined_lc)), bad_subclusters))



message("Cells before filtering: ", ncol(combined_lc)) # 5151
message("Cells after QC filtering: ", ncol(combined_lc_filtered)) # 4868




DimPlot(combined_lc_filtered,
        reduction = "umap", group.by  = "seurat_clusters",
        label = TRUE, repel = TRUE, pt.size = 0.4)


FeaturePlot(combined_lc_filtered,
            features  = c("nCount_RNA", "nFeature_RNA", "doublet_score"),
            reduction = "umap", pt.size   = 0.3,ncol      = 3)





saveRDS(combined_lc_filtered, output_file)

# 
# 

# 
# DimPlot(combined_lc_filtered, reduction = "umap",
#         group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
#   ggtitle("LC subclusters (filtered)")
# 
# # Relabel subclusters to 1..k (continuous)
# old_clusters <- sort(unique(Idents(combined_lc_filtered)))
# new_clusters <- seq_along(old_clusters)
# cluster_map <- setNames(as.character(new_clusters), old_clusters)
# 
# Idents(combined_lc_filtered) <- plyr::mapvalues(
#   Idents(combined_lc_filtered),
#   from = old_clusters,
#   to   = new_clusters
# )
# 
# # Store a dedicated column for subcluster ID
# combined_lc_filtered$subcluster <- Idents(combined_lc_filtered)
# 
# DimPlot(combined_lc_filtered, reduction = "umap",
#         group.by = "subcluster", label = TRUE, repel = TRUE) +
#   ggtitle("LC subclusters (renumbered)")

