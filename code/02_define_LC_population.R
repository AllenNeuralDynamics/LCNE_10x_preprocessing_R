# Select LC clusters from integrated object and save subset
source("./func/cluster_functions.R")

input_file  <- "../scratch/LC_clustered_harmony.rds"  # 231k cells
output_file <- "../scratch/LC_cluster_interest.rds"

# LC defined based on marker expression (see 01_integration_and_LC_selection.Rmd)
lc_markers <- c("Dbh", "Slc6a2", "Th", "Slc18a2", "Tacr3")



library(Seurat)

combined <- readRDS(input_file)
Idents(combined) <- "seurat_clusters"


# here we add the filteirng for this cluster 

lc_candidate_clusters <- get_candidate_clusters(combined, lc_markers[1:4])
combined_lc <- subset(combined, idents = lc_candidate_clusters) 
dim(combined_lc)  # 4984
dim(combined)  ## 231500
saveRDS(combined_lc, file = output_file)
