input_file = "../scratch/LC_cluster_interest.rds"  # 5151 cells to start with 
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
set.seed(42) 
combined_lc <- RunPCA(combined_lc, npcs = 30)
combined_lc <- RunUMAP(combined_lc, dims = 1:10, return.model = TRUE)

combined_lc <- FindNeighbors(combined_lc, dims = 1:10)
combined_lc <- FindClusters(combined_lc, resolution = 0.4,random.seed = 42)



# QC filtering steps

combined_lc$nCount_RNA   <- as.numeric(combined_lc$nCount_RNA)
combined_lc$nFeature_RNA <- as.numeric(combined_lc$nFeature_RNA)
combined_lc$percent.mt   <- as.numeric(combined_lc$percent.mt)

qc_df <- SeuratObject::FetchData(
  combined_lc,
  vars = c("nCount_RNA", "nFeature_RNA", "percent.mt", "seurat_clusters"),
  layer = "data"  # or whatever layer your QC metrics live in
  )
qc_df$seurat_clusters <- as.factor(qc_df$seurat_clusters)


combined_lc$nCount_RNA   <- as.numeric(combined_lc$nCount_RNA)
combined_lc$nFeature_RNA <- as.numeric(combined_lc$nFeature_RNA)
combined_lc$percent.mt   <- as.numeric(combined_lc$percent.mt)


# new: adding QC visualizations and filteirng 

qc_features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
qc_features <- qc_features[qc_features %in% colnames(combined_lc@meta.data)]
if ("doublet_score" %in% colnames(combined_lc@meta.data)) {
  qc_features <- c(qc_features, "doublet_score")
}
sapply(qc_features, function(f) class(combined_lc@meta.data[[f]]))
to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
for (f in qc_features) {
  combined_lc[[f]] <- to_num(combined_lc@meta.data[[f]])}

qc_df <- FetchData(combined_lc, vars = c("seurat_clusters", qc_features))
# qc_means <- aggregate(
#   qc_df[, qc_features],
#   by = list(cluster = qc_df$seurat_clusters),
#   FUN = mean,
#   na.rm = TRUE)
qc_medians <- aggregate(
  qc_df[, qc_features],
  by = list(cluster = qc_df$seurat_clusters),
  FUN = median,
  na.rm = TRUE)


# filter: the lowest detected RNA, and highest double score group are removed. 
qc_medians$cluster <- as.integer(as.character(qc_medians$cluster))
lowest_nCount  <- qc_medians$cluster[which.min(qc_medians$nCount_RNA)]
lowest_nFeature <- qc_medians$cluster[which.min(qc_medians$nFeature_RNA)]
highest_doublet <- qc_medians$cluster[which.max(qc_medians$doublet_score)]

rules <- c(
  nCount_RNA = "min",
  nFeature_RNA = "min",
  doublet_score = "max")
selected_clusters <- sapply(names(rules), function(col) {
  if (rules[col] == "min") {
    qc_medians$cluster[which.min(qc_medians[[col]])]
  } else {
    qc_medians$cluster[which.max(qc_medians[[col]])]}
})

selected_clusters

bad_subclusters <- as.character(unique(c(
  lowest_nCount,
  lowest_nFeature,
  highest_doublet
)))

bad_subclusters  # should be 5 and 7  


#####

combined_lc$cluster_badcluster_vs_other <- ifelse(
  combined_lc$seurat_clusters %in% bad_subclusters,
  "badcluster_candidates",
  "other"
)

combined_lc_filtered <- subset(combined_lc, idents = setdiff(levels(Idents(combined_lc)), bad_subclusters))


message("Cells before filtering: ", ncol(combined_lc)) # 4984
message("Cells after QC filtering: ", ncol(combined_lc_filtered)) # 4768

saveRDS(combined_lc_filtered, output_file)


