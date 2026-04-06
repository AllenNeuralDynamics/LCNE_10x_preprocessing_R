# Select LC clusters from integrated object and save subset

input_file  <- "../scratch/LC_clustered_harmony.rds"
output_file <- "../scratch/LC_cluster_interest.rds"

# LC defined based on marker expression (see 01_integration_and_LC_selection.Rmd)
lc_clusters_final <- c(9)   # 5151 samples will be selected

library(Seurat)

combined <- readRDS(input_file)
Idents(combined) <- "seurat_clusters"
combined_lc <- subset(combined, idents = lc_clusters_final) # 5151 cells -> 4984 (in mt version)
saveRDS(combined_lc, file = output_file)





