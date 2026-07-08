get_candidate_clusters <- function(obj, markers, thresh = 1, group.by = "seurat_clusters") {
  avg_mat <- AverageExpression(obj, features = markers, group.by = group.by,
                               slot = "data", verbose = FALSE)$RNA
  keep <- colnames(avg_mat)[colSums(avg_mat > thresh) == length(markers)]
  clusters <- as.integer(sub("^g", "", keep))
  
  n_cells <- table(obj@meta.data[[group.by]])[keep]
  cat(sprintf("Selected cluster %s: %d cells\n", clusters, n_cells), sep = "")
  
  clusters
}