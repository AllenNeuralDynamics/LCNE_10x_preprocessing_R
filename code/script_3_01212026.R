rm(list = ls())
gc()  # Garbage collection at start

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)

combined_lc <- readRDS("../scratch/LC_cluster_interest.rds")

Idents(combined_lc) <- "seurat_clusters"  # or another column if you set one

# Basic preprocessing within LC population
combined_lc <- NormalizeData(combined_lc)
combined_lc <- FindVariableFeatures(combined_lc, nfeatures = 2000)
combined_lc <- ScaleData(combined_lc)

# Dimensional reduction
combined_lc <- RunPCA(combined_lc, npcs = 30)
combined_lc <- RunUMAP(combined_lc, dims = 1:10, return.model = TRUE)

# Graph + clustering
combined_lc <- FindNeighbors(combined_lc, dims = 1:10)
combined_lc <- FindClusters(combined_lc, resolution = 0.4)

# Inspect subclusters
DimPlot(combined_lc, reduction = "umap", group.by = "seurat_clusters",
        label = TRUE, repel = TRUE, pt.size = 0.5) +
  ggtitle("LC population: subclusters")

# QC across subclusters
combined_lc$nCount_RNA   <- as.numeric(combined_lc$nCount_RNA)
combined_lc$nFeature_RNA <- as.numeric(combined_lc$nFeature_RNA)
combined_lc$percent.mt   <- as.numeric(combined_lc$percent.mt)
# Fetch data using the new API (no slot argument)
qc_df <- SeuratObject::FetchData(
  combined_lc,
  vars = c("nCount_RNA", "nFeature_RNA", "percent.mt", "seurat_clusters"),
  layer = "data"  # or whatever layer your QC metrics live in
)

qc_df$seurat_clusters <- as.factor(qc_df$seurat_clusters)

p1 <- ggplot(qc_df, aes(x = seurat_clusters, y = nCount_RNA)) +
  geom_violin() +
  theme_bw() +
  ggtitle("nCount_RNA")

p2 <- ggplot(qc_df, aes(x = seurat_clusters, y = nFeature_RNA)) +
  geom_violin() +
  theme_bw() +
  ggtitle("nFeature_RNA")

p3 <- ggplot(qc_df, aes(x = seurat_clusters, y = percent.mt)) +
  geom_violin() +
  theme_bw() +
  ggtitle("percent.mt")

p1
p2
p3
VlnPlot(combined_lc, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
group.by = "seurat_clusters", pt.size = 0.2, ncol = 3)










#--------------------------------------------
# Filter out low-quality subclusters and relabel
#--------------------------------------------

### -- new: adding analysis on the marker genes to check if we should remove 5 and 7 
Idents(combined_lc) <- "seurat_clusters"
lc_markers <- c("Dbh", "Slc6a2", "Th", "Slc18a2", "Tacr3")

FeaturePlot(
  combined_lc,
  features  = lc_markers,
  reduction = "umap",
  pt.size   = 0.3,
  ncol      = 3
)

# Also see the subclusters overlaid
DimPlot(
  combined_lc,
  reduction = "umap",
  group.by  = "seurat_clusters",
  label     = TRUE, repel = TRUE, pt.size = 0.4
)


# violin 
combined_lc$nCount_RNA   <- as.numeric(combined_lc$nCount_RNA)
combined_lc$nFeature_RNA <- as.numeric(combined_lc$nFeature_RNA)
combined_lc$percent.mt   <- as.numeric(combined_lc$percent.mt)

clusters_to_check <- c("0","2", "5", "7")

VlnPlot(
  subset(combined_lc, idents = clusters_to_check),
  features = lc_markers,
  group.by = "seurat_clusters",
  pt.size  = 0.1,
  ncol     = 3
)


## QC for these clusters
VlnPlot(
  subset(combined_lc, idents = clusters_to_check),
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  group.by = "seurat_clusters",
  pt.size  = 0.1,
  ncol     = 3
)

combined_lc$cluster_5_7_vs_other <- ifelse(
  combined_lc$seurat_clusters %in% c("5", "7"),
  "5_7_candidate",
  "other"
)

DotPlot(
  combined_lc,
  features = lc_markers,
  group.by = "cluster_5_7_vs_other"
) + RotatedAxis()




# remove 5 and 7 

bad_subclusters <- c("5","7")
combined_lc_filtered <- subset(combined_lc,
  idents = setdiff(levels(Idents(combined_lc)), bad_subclusters))
dim(combined_lc_filtered) # 4886 left remoivng 5, removing both 5 and 7 gives 4868 cells in total 


DimPlot(combined_lc_filtered, reduction = "umap",
        group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
  ggtitle("LC subclusters (filtered)")

# Relabel subclusters to 1..k (continuous)
old_clusters <- sort(unique(Idents(combined_lc_filtered)))
new_clusters <- seq_along(old_clusters)
cluster_map <- setNames(as.character(new_clusters), old_clusters)

Idents(combined_lc_filtered) <- plyr::mapvalues(
  Idents(combined_lc_filtered),
  from = old_clusters,
  to   = new_clusters
)

# Store a dedicated column for subcluster ID
combined_lc_filtered$subcluster <- Idents(combined_lc_filtered)

DimPlot(combined_lc_filtered, reduction = "umap",
        group.by = "subcluster", label = TRUE, repel = TRUE) +
  ggtitle("LC subclusters (renumbered)")

# Save for downstream work
dim(combined_lc_filtered)  #29617  4868
saveRDS(combined_lc_filtered, "../results/LC_subclusters_filtered.rds")
### the above corresponds to 'cluster_9_subclusters_filtered2.rds' from the original script!
