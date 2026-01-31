# Load required libraries

rm(list=ls())
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
# setwd('~/syoh/lc_seq_data/snrna_seq/')

file3= '~/syoh/lc_seq_data/snrna_seq/10xV4_Neuromodulatory_Noradrenergic_complete_RTX-4134_mat_241111.rda'
file1= '~/syoh/lc_seq_data/snrna_seq/10xV4_Neuromodulatory_Noradrenergic_RTX-4123_CR8_mat_240904.rda'
file2=  '~/syoh/lc_seq_data/snrna_seq/10xV4_Neuromodulatory_Noradrenergic_RTX-4132-4133_CR8_mat_241017.rda'

metafile3= '~/syoh/lc_seq_data/snrna_seq/10xV4_Neuromodulatory_Noradrenergic_complete_RTX-4134_CR8_samp.dat_241111.rda'
metafile2= '~/syoh/lc_seq_data/snrna_seq/10xV4_Neuromodulatory_Noradrenergic_RTX-4132-4133_CR8_samp.dat_241017.rda'
metafile1= '~/syoh/lc_seq_data/snrna_seq/10xV4_Neuromodulatory_Noradrenergic_RTX-4123_CR8_samp.dat_240904.rda'



# Processing each sample
#sample1 <- load_and_create_seurat(file1, metafile1, "Sample1", doublet_scores_col = "doublet_score")
#sample2 <- load_and_create_seurat(file2, metafile2, "Sample2", doublet_scores_col = "doublet_score")
sample3 <- load_and_create_seurat(file3, metafile3, "Sample3", doublet_scores_col = "doublet_score")


#saveRDS(sample1, "LC_sample1_filtered.rds")
#saveRDS(sample2, "LC_sample2_filtered.rds")

#it turns out that the third file has everything.

saveRDS(sample3, "LC_sample3_filtered.rds")

# Integrating the samples
sample_list <- list(Sample3 = sample3)
integrated_obj <- integrate_samples_hierarchical(sample_list, file_col = "batch_vendor_name", port_well_col = "rna_amplification")

# Visualizing integration results
integration_plots <- plot_integration_results(integrated_obj)
print(integration_plots$batch_effects)
print(integration_plots$clusters)


#-------------------------------------------------------------------------------
#functions
#-------------------------------------------------------------------------------
# Load required libraries
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)

load_and_create_seurat <- function(data_file, metadata_file, sample_name, 
                                   min_genes = 2000, max_genes = 15000, 
                                   max_mt_pct = 4, max_rb_pct = 4, 
                                   doublet_scores_col = 'doublet_score') {
  #data_file = file3
  #metadata_file = metafile3
  #sample_name = 'Sample3'
  
  # Load the data and metadata
  load(data_file)  # Assumes this loads a matrix of raw counts
  load(metadata_file)  # Assumes this loads a metadata data.frame
  
  # Filter metadata for rows under 'studies' belonging to 'Neuromodulatory_Noradrenergic'
  if (!"studies" %in% colnames(samp.dat)) {
    stop("The metadata file does not contain a 'studies' column.")
  }
  
  samp.dat <- samp.dat[samp.dat$studies == 'Neuromodulatory_Noradrenergic', ]
  
  # Ensure that mat and samp.dat dimensions match after filtering
  mat <- mat[, colnames(mat) %in% samp.dat$sample_id]
  if (ncol(mat) != nrow(samp.dat)) {
    stop("Counts matrix and metadata dimensions do not match after filtering!")
  }
  
  # Create Seurat object
  seurat_object <- CreateSeuratObject(
    counts = mat, 
    meta.data = samp.dat, 
    project = sample_name, 
    min.cells = 3, 
    min.features = 200
  )
  
  # Add metadata for experiment and platform
  seurat_object <- AddMetaData(seurat_object, metadata = sample_name, col.name = 'experiment')
  seurat_object <- AddMetaData(seurat_object, metadata = '10xV4', col.name = 'platform')
  
  # Calculate percentage of mitochondrial and ribosomal genes
  seurat_object$percent.mt <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  seurat_object$percent.rb <- PercentageFeatureSet(seurat_object, pattern = "^RP[SL]")
  
  # Add doublet scores if provided in the metadata
  if (!is.null(doublet_scores_col)) {
    if (doublet_scores_col %in% colnames(samp.dat)) {
      seurat_object$doublet_score <- samp.dat[[doublet_scores_col]]
    } else {
      warning(paste0("Doublet score column '", doublet_scores_col, "' not found in metadata. Skipping doublet score addition."))
    }
  }
  
  # Filter cells with QC metrics
  if (!is.null(doublet_scores_col) && doublet_scores_col %in% colnames(samp.dat)) {
    seurat_object <- subset(seurat_object, 
                            subset = nFeature_RNA > min_genes &
                              nFeature_RNA < max_genes &
                              percent.mt < max_mt_pct &
                              percent.rb < max_rb_pct &
                              doublet_score < 0.4)
  } else {
    seurat_object <- subset(seurat_object, 
                            subset = nFeature_RNA > min_genes &
                              nFeature_RNA < max_genes &
                              percent.mt < max_mt_pct &
                              percent.rb < max_rb_pct)
    warning("No doublet scores provided. Filtering only on count and QC metrics.")
  }
  
  # Normalize data
  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Find variable features
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 1500, assay = "RNA")
  
  rm(mat, samp.dat)
  
  return(seurat_object)
}




# Function to integrate samples using hierarchical batch correction
integrate_samples_hierarchical <- function(sample_list, file_col = "batch_vendor_name", port_well_col = "rna_amplification") {
  # Merge all samples
  #combined <- merge(x = sample_list[[1]], 
  #                  y = sample_list[2:length(sample_list)],
  #                  add.cell.ids = names(sample_list))
  # Replace all NA values in metadata with 'not_available'
  #load("LC_sample3_filtered.rds")
  metadata <- seurat_object@meta.data
  
  # Loop through each column in the metadata
  for (col in colnames(metadata)) {
    if (any(is.na(metadata[[col]]))) {
      metadata[[col]][is.na(metadata[[col]])] <- "not_available"
      message(paste("Replaced NA values in column:", col))
    } else {
      message(paste("No NA values in column:", col))
    }
  }
  
  # Update the Seurat object with modified metadata
  seurat_object@meta.data <- metadata
  
  # Verify that there are no NA values remaining
  if (any(is.na(seurat_object@meta.data))) {
    message("There are still NA values in the metadata.")
  } else {
    message("All NA values have been replaced with 'not_available'.")
  }
  
  # Scale data
  combined <- ScaleData(seurat_object, assay = "RNA")
  
  # Run PCA
  combined <- RunPCA(combined, assay = "RNA", npcs = 30)
  
  # Perform first level batch correction (e.g., by `batch_vendor_name`, RTX number)
  
  print("Running first-level Harmony batch correction (by RTX batch)...")
  combined <- RunHarmony(
    object = combined,         # Specify the Seurat object
    group.by.vars = file_col,  # The column containing batch/group information
    reduction = "pca",         # Base reduction to use
    assay.use = "RNA",         # Assay to use
    project.dim = FALSE,       # Whether to project onto dimensions
    verbose = TRUE             # Print progress
  )
  
  
  # Load Harmony library
  library(harmony)
  
  # Verify the batch column in metadata
  table(combined$batch)
  
  
  
  #-----------
  # Perform Harmony integration
  # Load Harmony library
  library(harmony)
  
  # Ensure the batch metadata is properly set
  table(combined$batch)
  combined$batch <- factor(combined$batch)
  
  # Ensure batch column exists and is properly set
  if (!"batch" %in% colnames(combined@meta.data)) {
    stop("The metadata column 'batch' does not exist in the Seurat object.")
  }
  
  
  # Step 1: Run Harmony for 'batch' correction
  combined <- RunHarmony(
    object = combined,
    reduction = "pca",
    group.by.vars = "batch",
    dims = 1:30
  )
  
  # Save Harmony results as the PCA embeddings
  combined@reductions$harmony_batch <- combined@reductions$harmony
  
  # Step 2: Run Harmony for 'rna_amplification' correction using batch-corrected embeddings
  # Set batch Harmony embeddings as the active reduction
  DefaultAssay(combined) <- "RNA"
  
  # Replace the PCA embeddings with Harmony batch embeddings for the next step
  combined <- RunHarmony(
    object = combined,
    group.by.vars = "rna_amplification",  # Column in metadata for RNA amplification information
    dims = 1:30                          # Number of dimensions to use
  )
  
  # Rename the Harmony reduction to reflect RNA amplification correction
  combined@reductions$harmony_final <- combined@reductions$harmony
  # Verify the reduction is stored
  print("Final Harmony reduction stored as 'harmony_final'.")
  
  # Add Harmony embeddings to the Seurat object
  combined <- RunUMAP(
    object = combined,
    reduction = "harmony",
    dims = 1:30
  )
  
  # Visualize the results
  DimPlot(combined, reduction = "umap", group.by = "batch",pt.size = 0.5)
  FeaturePlot(combined, features = "Dbh", reduction = "umap", pt.size = 0.5)
  FeaturePlot(combined, features = "Slc6a2", reduction = "umap", pt.size = 0.5)
  FeaturePlot(combined, features = "Th", reduction = "umap", pt.size = 0.5)
  FeaturePlot(combined, features = "Slc18a2", reduction = "umap", pt.size = 0.5)
  FeaturePlot(combined, features = "Tacr3", reduction = "umap", pt.size = 0.5)
  
  FeaturePlot(combined, features = "doublet_score", reduction = "umap", pt.size = 0.5, cols = c("lightgrey", "red"))
  
  
  # Find neighbors using final Harmony embeddings
  combined <- FindNeighbors(combined, 
                            reduction = "harmony",
                            dims = 1:30)
  
  # Find clusters
  combined <- FindClusters(combined, resolution = 0.5, random.seed = 42)
  
  return(combined)
}

# Function to visualize integration results
plot_integration_results <- function(seurat_obj) {
  p1 <- DimPlot(combined, reduction = "pca", group.by = "orig.ident", pt.size = 0.1) +
    ggtitle("Before Harmony Integration")
  p2 <- DimPlot(combined, reduction = "harmony", group.by = "orig.ident", pt.size = 0.1) +
    ggtitle("After Harmony Integration")
  p3 <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 0.1) +
    ggtitle("Clusters")
  return(list(batch_effects = p1 + p2, clusters = p3))
}

#-------------------------------
#Focus on cluser 9 and 24
#-------------------------------
# Subset the Seurat object for clusters 9 and 24
subset_clusters <- subset(combined, idents = c(9, 24))

# Visualize read counts for clusters 9 and 24
VlnPlot(subset_clusters, features = "nCount_RNA", pt.size = 0.5) +
  theme_minimal() + 
  ggtitle("Read Counts for Clusters 9 and 24")

saveRDS(combined, "LC_clustered_harmony.rds")
#looks like cluster 24 is safely junk
combined <- readRDS("LC_clustered_harmony.rds")
#--------------------------------
#Focus on cluster 9
#------------------------------
install.packages('devtools')
#the following did not work in HPC
devtools::install_github('immunogenomics/presto')

cat("Number of cells in Cluster 9:", ncol(cluster_9), "\n")
summary(cluster_9@meta.data)
VlnPlot(combined, features = c("nCount_RNA", "percent.mt"), group.by = "seurat_clusters", idents = 9)
cluster_9_markers <- FindMarkers(combined, ident.1 = 9)
head(cluster_9_markers)


#-------------------------------
#analysis of cluster 9 starts here
#-------------------------------

# Subset cluster 9
cluster_9 <- subset(combined, idents = 9)
#saveRDS(cluster_9, "cluster_9.rds")

# Normalize and scale data
cluster_9 <- NormalizeData(cluster_9)
cluster_9 <- FindVariableFeatures(cluster_9) #default is 2000, maybe 1500 is better
cluster_9 <- ScaleData(cluster_9)

# Perform dimensional reduction
cluster_9 <- RunPCA(cluster_9, npcs = 30)
cluster_9 <- RunUMAP(cluster_9, dims = 1:10, return.model = T)  # Adjust dimensions as needed

# Find neighbors and clusters
cluster_9 <- FindNeighbors(cluster_9, dims = 1:10)
cluster_9 <- FindClusters(cluster_9, resolution = 0.4)  # Adjust resolution as needed

# Visualize subclusters
DimPlot(cluster_9, reduction = "umap", group.by = "seurat_clusters", label = F, pt.size = 0.5) + 
  ggtitle("Subclusters of Cluster 9")

# Analyze subclusters - Find marker genes
subcluster_markers <- FindAllMarkers(cluster_9)
print("Top markers for each subcluster:")
head(subcluster_markers)

# Get cell statistics for subclusters
subcluster_stats <- table(Idents(cluster_9))
print("Number of cells in each subcluster:")
print(subcluster_stats)

# (Optional) Add subcluster information back to the original Seurat object
combined$subcluster_9 <- "not_cluster_9"
combined$subcluster_9[Cells(cluster_9)] <- Idents(cluster_9)

# (Optional) Visualize subclusters on the original UMAP
DimPlot(combined, reduction = "umap", group.by = "subcluster_9", label = TRUE) +
  ggtitle("Subcluster 9 Integration")


write.csv(subcluster_markers, "cluster9_markers.csv")

saveRDS(cluster_9, "cluster_9_subclusters.rds")
saveRDS(combined, "LC_harmony_cluster_with_cluster_9_subclusters.rds")

#starting in the middle
cluster_9 <- readRDS("cluster_9_subclusters.rds")
subcluster_markers <- read.csv("cluster9_markers.csv", stringsAsFactors = F)
cluster_9 <- RunUMAP(cluster_9, dims = 1:10, return.model = T)  # Adjust dimensions as needed

# Verify the filtering
DimPlot(cluster_9, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Filtered Subclusters of Cluster 9 (Excluding Clusters 5 and 7)")

FeaturePlot(cluster_9, features = "Pald1", reduction = "umap", pt.size = 0.5)

#------------------------------------------------------------
#Plotting
#------------------------------------------------------------
#filter bad clusters:
#first plot histograms of read counts among the subclusters
#total gene counts as well
# Visualize read counts for clusters 9 and 24
VlnPlot(cluster_9, features = "nFeature_RNA", pt.size = 0.5) +
  theme_minimal() + 
  ggtitle("Read Counts for Clusters 9")

#cluster 5 and 7 are weird - throw away


# Exclude clusters 5 and 7
cluster_9_filtered <- subset(cluster_9, idents = setdiff(unique(Idents(cluster_9)), c(5, 7)))

# Verify the filtering
DimPlot(cluster_9_filtered, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Filtered Subclusters of Cluster 9 (Excluding Clusters 5 and 7)")

## Rename clusters to start from 1 and be continuous
old_clusters <- unique(Idents(cluster_9_filtered))
new_clusters <- seq_along(old_clusters)  # Create a continuous sequence
cluster_mapping <- setNames(new_clusters, old_clusters)  # Map old to new

# Apply the new cluster names
Idents(cluster_9_filtered) <- factor(Idents(cluster_9_filtered), levels = old_clusters)
Idents(cluster_9_filtered) <- plyr::mapvalues(Idents(cluster_9_filtered), from = old_clusters, to = new_clusters)

# Update the seurat_clusters column in the metadata
cluster_9_filtered$seurat_clusters <- Idents(cluster_9_filtered)

# Ensure the identities are correctly set
Idents(cluster_9_filtered) <- cluster_9_filtered$seurat_clusters

# Plot with the updated cluster names
DimPlot(cluster_9_filtered, reduction = "umap", group.by = "seurat_clusters", label = F, pt.size = 0.6) +
  ggtitle("Renamed Subclusters of Cluster 9")

saveRDS(cluster_9_filtered, "cluster_9_subclusters_filtered2.rds")

#Starting from the middle
cluster_9_filtered <-readRDS("cluster_9_subclusters_filtered2.rds")

#-------------------------------------------------------
#Plot like in python
#-------------------------------------------------------
library(ggplot2)
library(RColorBrewer)

# Extract UMAP embeddings and cluster labels
umap_data <- as.data.frame(Embeddings(cluster_9_filtered, reduction = "umap"))
umap_data$cluster <- Idents(cluster_9_filtered)

# Ensure consistent column names
colnames(umap_data)[1:2] <- c("umap_1", "umap_2")

# Choose a palette with enough distinct colors
n_clusters <- length(unique(umap_data$cluster))
palette_colors <- brewer.pal(min(n_clusters, 8), "Dark2")
if (n_clusters > 8) palette_colors <- colorRampPalette(palette_colors)(n_clusters)

# Plot with enhanced aesthetics
ggplot(umap_data, aes(x = umap_1, y = umap_2, color = factor(cluster))) +
  geom_point(size = 1.5, alpha = 0.4, shape = 16) +  # Slightly larger, opaque points
  scale_color_manual(values = palette_colors) +      # Improved color palette
  labs(
    title = "UMAP of Cluster 9 Subpopulations",
    subtitle = "Colored by subcluster identity",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2",
    color = "Subcluster"
  ) +
  theme_bw(base_size = 14) +  # Clean, professional theme
  guides(color = guide_legend(override.aes = list(size = 4))) +  # Enlarged legend points
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Centered title
    plot.subtitle = element_text(hjust = 0.5, size = 12),             # Subtitle
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.position = "right",                                       # Legend on the side
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90"),               # Subtle major gridlines
    panel.grid.minor = element_blank(),                              # No minor gridlines
    panel.border = element_rect(color = "black", fill = NA, size = 0.8)  # Add a border
  )


#-----------------------------------------------------------------
#Changing cluster number so that it can go from top to bottom
#-----------------------------------------------------------------
library(ggplot2)
library(RColorBrewer)

# Extract UMAP embeddings and cluster labels
umap_data <- as.data.frame(Embeddings(cluster_9_filtered, reduction = "umap"))
umap_data$cluster <- Idents(cluster_9_filtered)

# Ensure consistent column names
colnames(umap_data)[1:2] <- c("umap_1", "umap_2")

# Create a mapping of old to new cluster names
cluster_mapping <- c("6" = "1", "2" = "2", "4" = "3", "3" = "4", "5" = "5", "1" = "6")

# Apply the mapping to rename clusters
umap_data$cluster <- factor(cluster_mapping[as.character(umap_data$cluster)], levels = c("1", "2", "3", "4", "5", "6"))

# Define fixed colors based on the old cluster names
original_colors <- brewer.pal(8, "Dark2")
fixed_color_mapping <- c(
  "1" = original_colors[6],  # Old cluster 6 -> New cluster 1
  "2" = original_colors[2],  # Old cluster 2 -> New cluster 2
  "3" = original_colors[4],  # Old cluster 4 -> New cluster 3
  "4" = original_colors[3],  # Old cluster 3 -> New cluster 4
  "5" = original_colors[5],  # Old cluster 5 -> New cluster 5
  "6" = original_colors[1]   # Old cluster 1 -> New cluster 6
)

# Plot with renamed clusters and preserved colors
ggplot(umap_data, aes(x = umap_1, y = umap_2, color = factor(cluster))) +
  geom_point(size = 1.5, alpha = 0.4, shape = 16) +  # Main plot points
  scale_color_manual(values = fixed_color_mapping) +  # Fixed color mapping
  labs(
    title = "UMAP of Cluster 9 Subpopulations (Renamed)",
    subtitle = "Colored by subcluster identity",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2",
    color = "LC-NE cluster"
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) +  # Enlarged legend points
  theme_bw(base_size = 14) +  # Clean theme
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Centered title
    plot.subtitle = element_text(hjust = 0.5, size = 12),             # Subtitle
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.position = "right",                                       # Legend on the side
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90"),               # Subtle gridlines
    panel.grid.minor = element_blank(),                              # No minor gridlines
    panel.border = element_rect(color = "black", fill = NA, size = 0.8)  # Border
  )




#----------------------------------------------------------------------------------
#Draw markers
#----------------------------------------------------------------------------------
# Load necessary libraries
library(dplyr)

# Load the original marker gene file
cluster9_markers <- read.csv("cluster9_markers.csv")

# Step 1: Remove rows corresponding to clusters 5 and 7
filtered_markers <- cluster9_markers %>%
  filter(!cluster %in% c(5, 7))

# Step 2: Map remaining clusters to sequential names (low-to-high mapping)
low_to_high_mapping <- c("0" = "1", "1" = "2", "2" = "3", "3" = "4", "4" = "5", "6" = "6")
filtered_markers <- filtered_markers %>%
  mutate(IntermediateCluster = low_to_high_mapping[as.character(cluster)])

# Step 3: Apply the final mapping
final_mapping <- c("6" = "1", "2" = "2", "4" = "3", "3" = "4", "5" = "5", "1" = "6")
filtered_markers <- filtered_markers %>%
  mutate(UpdatedCluster = final_mapping[IntermediateCluster])

# Step 4: Save the updated marker gene file
write.csv(filtered_markers, "cluster9_markers_updated.csv", row.names = FALSE)

cat("Updated marker gene file saved as 'cluster9_markers_updated.csv'\n")


#------------
#plot genes
#------------
#need to plot canonical marker genes: Th, Dbh, Trhr, Pdyn, Gpr101, Tacr3

genes_to_plot = c("Th", "Dbh", "Trhr", "Pdyn", "Gpr101", "Tacr3")

#check out what the expression unit is
str(GetAssayData(cluster_9_filtered, assay = "RNA", layer = "data"))
hist(summary(GetAssayData(cluster_9_filtered, assay = "RNA", layer = "data"))$x)
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Load updated marker gene file
marker_genes <- read.csv("cluster9_markers_updated.csv")

# Specify the cluster to plot
cluster_to_plot <- "6"  # Replace with the desired cluster number as a string

# Identify the top 20 marker genes for the selected cluster based on fold change
top_genes <- marker_genes %>%
  filter(UpdatedCluster == cluster_to_plot) %>%
  arrange(desc(avg_log2FC)) %>%  # Arrange by descending fold change
  slice_head(n = 20) %>%  # Select top 20
  pull(gene)

# Load UMAP embeddings
umap_data <- as.data.frame(Embeddings(cluster_9_filtered, reduction = "umap"))
umap_data$Cell <- rownames(umap_data)  # Add cell IDs for merging

# Load normalized expression data
# Replace "RNA" with the assay name in your Seurat object if different
expr_data <- as.data.frame(GetAssayData(cluster_9_filtered, assay = "RNA", layer = "data"))
expr_data <- expr_data[top_genes, , drop = FALSE]  # Subset for top genes
expr_data <- t(expr_data)  # Transpose so genes are columns
expr_data <- as.data.frame(expr_data)
expr_data$Cell <- rownames(expr_data)  # Add cell IDs for merging

# Merge UMAP and expression data
plot_data <- umap_data %>%
  left_join(expr_data, by = "Cell") %>%
  pivot_longer(cols = all_of(top_genes), names_to = "Gene", values_to = "Expression")

# Sort data points so that cells with higher expression are plotted last
plot_data <- plot_data %>%
  arrange(Gene, Expression)

# Plot UMAP with expression levels
ggplot(plot_data, aes(x = umap_1, y = umap_2, color = Expression)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_gradient(low = "grey80", high = "red") +  # Grey for low, red for high expression
  facet_wrap(~ Gene, ncol = 5) +  # Create separate panels for each gene
  theme_minimal(base_size = 12) +
  labs(
    title = paste("Top 20 Marker Genes for Cluster", cluster_to_plot),
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2",
    color = "Expression Level"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    strip.text = element_text(face = "bold", size = 10),  # Gene labels
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

#looks like I have to carefully curate what makes sense to show - p-value has interesting genes for cluster 1 but not other clusters.

#---------------------------
#plot random genes
#---------------------------
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Define genes of interest
genes_of_interest <- c("Tacr3", "Pdyn", "Trhr", "Shox2", "Dbh", "Npy", "Scn9a", "Cacna1i", "Kcnq5", "Cacna1i", "Gpr101")

# Load UMAP embeddings
umap_data <- as.data.frame(Embeddings(cluster_9_filtered, reduction = "umap"))
umap_data$Cell <- rownames(umap_data)  # Add cell IDs for merging

# Load normalized expression data
# Replace "RNA" with the assay name in your Seurat object if different
expr_data <- as.data.frame(GetAssayData(cluster_9_filtered, assay = "RNA", layer = "data"))

# Ensure all genes of interest are in the dataset
genes_to_plot <- intersect(genes_of_interest, rownames(expr_data))

if (length(genes_to_plot) < length(genes_of_interest)) {
  warning("Some genes of interest were not found in the dataset.")
}

# Subset expression data for the selected genes
expr_data <- expr_data[genes_to_plot, , drop = FALSE]
expr_data <- t(expr_data)  # Transpose so genes are columns
expr_data <- as.data.frame(expr_data)
expr_data$Cell <- rownames(expr_data)  # Add cell IDs for merging

# Merge UMAP and expression data
plot_data <- umap_data %>%
  left_join(expr_data, by = "Cell") %>%
  pivot_longer(cols = all_of(genes_to_plot), names_to = "Gene", values_to = "Expression")

# Sort data points so that cells with higher expression are plotted last
plot_data <- plot_data %>%
  arrange(Gene, Expression)

# Plot UMAP with expression levels
ggplot(plot_data, aes(x = umap_1, y = umap_2, color = Expression)) +
  geom_point(size = 1.5, alpha = 0.4) +
  scale_color_gradient(low = "grey95", high = "red") +  # Grey for low, red for high expression
  facet_wrap(~ Gene, ncol = 3) +  # Create separate panels for each gene
  theme_minimal(base_size = 12) +
  labs(
    title = "Expression of Genes of Interest in UMAP Space",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2",
    color = "Expression Level"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    strip.text = element_text(face = "bold", size = 10),  # Gene labels
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )


#plot separately

# Create a directory to save the plots
output_dir <- "gene_expression_plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Loop over each gene and generate a separate plot
for (gene in genes_to_plot) {
  # Subset data for the current gene
  gene_data <- plot_data %>%
    filter(Gene == gene)
  
  # Generate the plot
  p <- ggplot(gene_data, aes(x = umap_1, y = umap_2, color = Expression)) +
    geom_point(size = 1.5, alpha = 0.4) +
    scale_color_gradient(low = "grey95", high = "red") +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("Expression of", gene, "in UMAP Space"),
      x = "UMAP Dimension 1",
      y = "UMAP Dimension 2",
      color = "Expression Level"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  # Save the plot as a PNG file
  output_file <- file.path(output_dir, paste0(gene, "_expression.png"))
  ggsave(output_file, plot = p, width = 5, height = 5, dpi = 300)
}

cat("Plots saved in the directory:", output_dir)


#----------------------------------------
#More QC plots
#----------------------------------------
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tibble)

# Extract UMAP embeddings
umap_data <- as.data.frame(Embeddings(cluster_9_filtered, reduction = "umap"))
umap_data$Cell <- rownames(umap_data)  # Add cell IDs for merging

# Extract the 'rna_amplification' column and attach rownames as 'Cell'
rna_amplification_data <- cluster_9_filtered@meta.data %>%
  rownames_to_column(var = "Cell") %>%  # Convert rownames to a column
  select(Cell, rna_amplification)  # Select only the relevant column

# Merge UMAP data with 'rna_amplification'
plot_data <- umap_data %>%
  left_join(rna_amplification_data, by = "Cell")


# Generate a larger color palette
# Use RColorBrewer's "Set2" for up to 8 colors, and scale to more if needed
palette_size <- length(unique(plot_data$rna_amplification))
custom_colors <- if (palette_size <= 8) {
  RColorBrewer::brewer.pal(n = palette_size, name = "Set2")
} else {
  grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(palette_size)
}

# Plot UMAP with 'rna_amplification' as color
ggplot(plot_data, aes(x = umap_1, y = umap_2, color = rna_amplification)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = custom_colors) +
  theme_minimal(base_size = 12) +
  labs(
    title = "UMAP Colored by RNA Amplification",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2",
    color = "RNA Amplification"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )



#color code by sex


# Load necessary libraries
library(ggplot2)
library(dplyr)

# Extract UMAP embeddings
umap_data <- as.data.frame(Embeddings(cluster_9_filtered, reduction = "umap"))
colnames(umap_data)[1:2] <- c("umap_1", "umap_2")  # Rename columns for consistency
umap_data$Cell <- rownames(umap_data)  # Add cell IDs for merging

# Extract the 'sex' column and attach rownames as 'Cell'
sex_data <- cluster_9_filtered@meta.data %>%
  rownames_to_column(var = "Cell") %>%  # Convert rownames to a column
  select(Cell, sex)  # Select only the relevant column

# Merge UMAP data with 'sex'
plot_data <- umap_data %>%
  left_join(sex_data, by = "Cell")

# Generate a custom color palette for 'sex'
palette_size <- length(unique(plot_data$sex))
custom_colors <- if (palette_size <= 8) {
  RColorBrewer::brewer.pal(n = palette_size, name = "Set2")
} else {
  grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(palette_size)
}

# Plot UMAP with 'sex' as color
ggplot(plot_data, aes(x = umap_1, y = umap_2, color = sex)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = custom_colors) +
  theme_minimal(base_size = 12) +
  labs(
    title = "UMAP Colored by Sex",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2",
    color = "Sex"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

#potentially there is a sexually dimorphic cluster!

#----------------------------------------------------------------------------------
#To do: analysis of cluster 4 vs 5 - sexually dimorphic cluster or an artefact?
#----------------------------------------------------------------------------------

#check out autosomal genes that are enriched in one of the clusters
#the expression levels should be significantly different - normally about 2 fold


# List of female enriched genes from Mulvey et al. Cell Reports 2018
genes <- c(
  "Fam171b", "Gabre", "Cd200", "Dnaja3", "Tcf12", "Nat14", "Tead1", "Hbs1l", "Fam120b",
  "Senp3", "Gga2", "Dgkb", "Vrk1", "NA", "Tacr3", "Ptrh2", "Sf3b2", "NA", "Srd5a3",
  "Dguok", "Gipc1", "Eif2s2", "Slc6a15", "Ispd", "Wdr92", "Snhg12", "Tm7sf3", "Zfp92",
  "Lin28b", "Pja1", "Atp6v1a", "Raph1", "Zyg11b", "Ubxn7", "Peg10", "Nt5dc2", "Gosr2",
  "Nkrf", "NA", "Sgip1", "Polm", "Pacrg", "Tbc1d19", "Golt1b", "Iqcb1", "2310035C23Rik",
  "Mrpl3", "Arg2", "Mrps11", "Nfu1", "Arfgap1", "Gga2", "Stx17", "Stx17", "Elavl4",
  "Tiam1", "H13", "Ap1g1", "Afg3l1", "Sec22a", "Ndufs3", "Adgrb3", "Dzip1l", "Ptger3",
  "Disp1", "Stk32c", "Tex9", "Ntng1", "Cd151", "Fpgs", "Tcea3", "Zfp386"
)

# Remove "NA" and create a dataframe
filtered_genes <- genes[genes != "NA"]
df <- data.frame(Gene = filtered_genes)

# Print the dataframe
print(df)
DimPlot(snrnaseq_obj)
#plot
FeaturePlot(snrnaseq_obj, features = genes, reduction = "umap", ncol = 5) + 
  ggtitle("UMAP Visualization of Selected Genes") +
  theme_minimal()

#does not look like these are any different

# Subset Seurat object to include only clusters 3 and 5
cluster_subset <- subset(snrnaseq_obj, idents = c(3, 5))

# Calculate marker genes distinguishing cluster 3 from cluster 5
markers <- FindMarkers(cluster_subset, ident.1 = 3, ident.2 = 5, test.use = "wilcox")

# View top markers
head(markers)

# Filter significant markers (adjust p-value < 0.05)
significant_markers <- markers[markers$p_val_adj < 0.05, ]

# View top significant markers
significant_markers <- significant_markers[order(significant_markers$avg_log2FC, decreasing = TRUE), ]
head(significant_markers)

# Select the top 10 marker genes based on adjusted p-value and log fold change
top_markers <- rownames(significant_markers)[11:20]
print(top_markers)

# Plot UMAP for the top 10 marker genes
FeaturePlot(snrnaseq_obj, features = top_markers, reduction = "umap", ncol = 5) + 
  theme_minimal()
#looks like there isn't a good case for sexual dimorphic cluster - the male vs female markers split down the middle. Cluster 5 is in female but is defined by lack of markers

#--------------------------------------------------
#Plot in PCA space
#--------------------------------------------------

library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Extract PCA embeddings
pca_data <- as.data.frame(Embeddings(cluster_9_filtered, reduction = "pca"))
colnames(pca_data)[1:3] <- c("PCA_1", "PCA_2", "PCA_3")  # Rename PCA columns
pca_data$Cell <- rownames(pca_data)  # Add cell IDs for merging

# Extract cluster assignments
cluster_data <- cluster_9_filtered@meta.data %>%
  rownames_to_column(var = "Cell") %>%  # Convert rownames to a column
  select(Cell, cluster = seurat_clusters)  # Replace 'seurat_clusters' if another column contains clusters

# Create a mapping of old to new cluster names (ensure it's consistent with UMAP)
cluster_mapping <- c("6" = "1", "2" = "2", "4" = "3", "3" = "4", "5" = "5", "1" = "6")

# Apply the mapping to rename clusters
cluster_data$cluster <- factor(cluster_mapping[as.character(cluster_data$cluster)], levels = c("1", "2", "3", "4", "5", "6"))

# Merge PCA data with cluster information
plot_data <- pca_data %>%
  left_join(cluster_data, by = "Cell")

# Use the fixed color mapping from your UMAP plot
original_colors <- brewer.pal(8, "Dark2")
fixed_color_mapping <- c(
  "1" = original_colors[6],  # Old cluster 6 -> New cluster 1
  "2" = original_colors[2],  # Old cluster 2 -> New cluster 2
  "3" = original_colors[4],  # Old cluster 4 -> New cluster 3
  "4" = original_colors[3],  # Old cluster 3 -> New cluster 4
  "5" = original_colors[5],  # Old cluster 5 -> New cluster 5
  "6" = original_colors[1]   # Old cluster 1 -> New cluster 6
)

# Plot PCA 1 vs PCA 2
ggplot(plot_data, aes(x = PCA_1, y = PCA_2, color = factor(cluster))) +
  geom_point(size = 1.5, alpha = 0.6, shape = 16) +
  scale_color_manual(values = fixed_color_mapping) +
  labs(
    title = "PCA of Cluster 9 Subpopulations (Renamed)",
    x = "PCA Dimension 1",
    y = "PCA Dimension 2",
    color = "LC-NE cluster"
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10)
  )

# Use the fixed color mapping from your UMAP plot
library(plotly)

# 3D Scatter Plot for PCA 1, 2, 3
plot_ly(
  data = plot_data,
  x = ~PCA_1, y = ~PCA_2, z = ~PCA_3,
  color = ~factor(cluster),
  colors = fixed_color_mapping,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3, opacity = 0.8)
) %>%
  layout(
    title = "3D PCA of Cluster 9 Subpopulations",
    scene = list(
      xaxis = list(title = "PCA 1"),
      yaxis = list(title = "PCA 2"),
      zaxis = list(title = "PCA 3")
    )
  )


#--------------------------
#find a polynomial fit
#--------------------------

library(ggplot2)
library(dplyr)

# Fit a polynomial regression model (quadratic: degree 2)
poly_fit <- lm(PCA_2 ~ poly(PCA_1, 2), data = plot_data)

# Generate predictions for the curve
curve_data <- data.frame(PCA_1 = seq(min(plot_data$PCA_1), max(plot_data$PCA_1), length.out = 100))
curve_data$PCA_2 <- predict(poly_fit, newdata = curve_data)

# Plot PCA points and fitted curve
ggplot(plot_data, aes(x = PCA_1, y = PCA_2, color = factor(cluster))) +
  geom_point(size = 1.5, alpha = 0.6, shape = 16) +
  scale_color_manual(values = fixed_color_mapping) +
  geom_line(data = curve_data, aes(x = PCA_1, y = PCA_2), color = "black", size = 1) +
  labs(
    title = "PCA of Cluster 9 Subpopulations with Fitted Curve",
    x = "PCA Dimension 1",
    y = "PCA Dimension 2",
    color = "Cluster"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10)
  )


# in 3d

library(plotly)
library(dplyr)
#install.packages("princurve")
library(princurve)


# Extract PCA coordinates
pca_coords <- plot_data %>%
  select(PCA_1, PCA_2, PCA_3) %>%
  as.matrix()

# Fit a principal curve through the PCA points
fit <- principal_curve(pca_coords)

# Extract the fitted curve and rename columns
fitted_curve <- as.data.frame(fit$s[order(fit$lambda), ])
colnames(fitted_curve) <- c("PCA_1", "PCA_2", "PCA_3")  # Rename columns for clarity

# Add fitted curve to plotly visualization
plot_ly() %>%
  # Scatter points for PCA
  add_trace(
    data = plot_data,
    x = ~PCA_1, y = ~PCA_2, z = ~PCA_3,
    type = "scatter3d",
    mode = "markers",
    color = ~factor(cluster),
    colors = fixed_color_mapping,
    marker = list(size = 3, opacity = 0.8)
  ) %>%
  # Fitted curve in 3D
  add_trace(
    data = fitted_curve,
    x = ~PCA_1, y = ~PCA_2, z = ~PCA_3,
    type = "scatter3d",
    mode = "lines",
    line = list(color = "black", width = 4),
    name = "Fitted Curve"
  ) %>%
  layout(
    title = "3D PCA with Fitted Curve",
    scene = list(
      xaxis = list(title = "PCA 1"),
      yaxis = list(title = "PCA 2"),
      zaxis = list(title = "PCA 3")
    )
  )


#unroll

library(plotly)
library(dplyr)
library(princurve)

# Extract PCA coordinates
pca_coords <- plot_data %>%
  select(PCA_1, PCA_2, PCA_3) %>%
  as.matrix()

# Fit a principal curve
fit <- principal_curve(pca_coords)

# Extract the arc length (lambda) and the fitted points on the curve
plot_data$arc_length <- fit$lambda  # Unrolled position (arc length along the curve)

# Calculate distances to the curve
fitted_coords <- as.data.frame(fit$s)
colnames(fitted_coords) <- c("PCA_1", "PCA_2", "PCA_3")  # Name columns
plot_data$distance_to_curve <- sqrt(
  (plot_data$PCA_1 - fitted_coords$PCA_1)^2 +
    (plot_data$PCA_2 - fitted_coords$PCA_2)^2 +
    (plot_data$PCA_3 - fitted_coords$PCA_3)^2
)

# Visualize the unrolled arc in 2D
ggplot(plot_data, aes(x = arc_length, y = distance_to_curve, color = factor(cluster))) +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_color_manual(values = fixed_color_mapping) +
  labs(
    title = "Unrolled PCA Arc",
    x = "Arc Length (Unrolled Position)",
    y = "Distance to Curve",
    color = "Cluster"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10)
  )

#PCA loading


pca_loadings <- Loadings(cluster_9_filtered, reduction = "pca")
str(pca_loadings)

library(ggplot2)
library(dplyr)
library(tidyr)

# Convert PCA loadings to a data frame
pca_loadings <- as.data.frame(Loadings(cluster_9_filtered, reduction = "pca"))

# Add gene names as a column
pca_loadings$Gene <- rownames(pca_loadings)

# Analyze top contributing genes for PCA 1, 2, and 3
top_genes <- pca_loadings %>%
  select(Gene, PC_1, PC_2, PC_3) %>%  # Select loadings for PC 1, 2, and 3
  pivot_longer(cols = starts_with("PC"), names_to = "Principal_Component", values_to = "Loading") %>%
  group_by(Principal_Component) %>%
  arrange(desc(abs(Loading)), .by_group = TRUE) %>%
  slice_head(n = 100)  # Select top 100 genes for each principal component

# Save results to CSV file
write.csv(top_genes, "top_100_genes_pca1_3.csv", row.names = FALSE)

# Generate separate plots for PCA 1, PCA 2, and PCA 3
pc1_plot <- ggplot(top_genes %>% filter(Principal_Component == "PC_1"),
                   aes(x = reorder(Gene, abs(Loading)), y = abs(Loading))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 100 Contributing Genes to PC 1",
    x = "Gene",
    y = "Absolute Loading"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.text.y = element_text(size = 8)
  )

pc2_plot <- ggplot(top_genes %>% filter(Principal_Component == "PC_2"),
                   aes(x = reorder(Gene, abs(Loading)), y = abs(Loading))) +
  geom_bar(stat = "identity", fill = "tomato") +
  coord_flip() +
  labs(
    title = "Top 100 Contributing Genes to PC 2",
    x = "Gene",
    y = "Absolute Loading"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.text.y = element_text(size = 8)
  )

pc3_plot <- ggplot(top_genes %>% filter(Principal_Component == "PC_3"),
                   aes(x = reorder(Gene, abs(Loading)), y = abs(Loading))) +
  geom_bar(stat = "identity", fill = "forestgreen") +
  coord_flip() +
  labs(
    title = "Top 100 Contributing Genes to PC 3",
    x = "Gene",
    y = "Absolute Loading"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.text.y = element_text(size = 8)
  )

# Display plots
print(pc1_plot)
print(pc2_plot)
print(pc3_plot)

# Optionally save the plots to files
ggsave("PC1_Top_100_Genes.png", plot = pc1_plot, width = 10, height = 8, dpi = 300)
ggsave("PC2_Top_100_Genes.png", plot = pc2_plot, width = 10, height = 8, dpi = 300)
ggsave("PC3_Top_100_Genes.png", plot = pc3_plot, width = 10, height = 8, dpi = 300)

#dot plot
library(ggplot2)
library(dplyr)
library(tidyr)

# Convert PCA loadings to a data frame
pca_loadings <- as.data.frame(Loadings(cluster_9_filtered, reduction = "pca"))

# Add gene names as a column
pca_loadings$Gene <- rownames(pca_loadings)

# Analyze top contributing genes for PCA 1, 2, and 3
top_genes <- pca_loadings %>%
  select(Gene, PC_1, PC_2, PC_3) %>%  # Select loadings for PC 1, 2, and 3
  mutate(
    Absolute_PC_1 = abs(PC_1),
    Absolute_PC_2 = abs(PC_2),
    Absolute_PC_3 = abs(PC_3)
  )

# Sort genes by PC1, then PC2, and finally PC3 (all high to low)
sorted_genes <- top_genes %>%
  arrange(desc(Absolute_PC_1), desc(Absolute_PC_2), desc(Absolute_PC_3)) %>%
  pull(Gene)

# Keep only the top 100 genes across all PCs
top_genes_filtered <- top_genes %>%
  filter(Gene %in% sorted_genes[1:100]) %>%
  pivot_longer(cols = starts_with("PC"), names_to = "Principal_Component", values_to = "Loading") %>%
  mutate(Absolute_Loading = abs(Loading))

# Create dot plot
ggplot(top_genes_filtered, aes(
  x = Principal_Component,
  y = factor(Gene, levels = sorted_genes[1:100]),  # Sort by hierarchical order
  size = Absolute_Loading
)) +
  geom_point(color = "steelblue", alpha = 0.7) +
  scale_size_continuous(range = c(1, 10), name = "Absolute Loading") +  # Adjust circle size
  labs(
    title = "Top 100 Contributing Genes to PC 1, 2, and 3 (Hierarchically Sorted)",
    x = "Principal Component",
    y = "Gene"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10)
  )


#heatmap

library(ggplot2)
library(dplyr)
library(tidyr)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
#install.packages("colorRamp2")
library(colorRamp2)
library(circlize)


# Convert PCA loadings to a data frame
pca_loadings <- as.data.frame(Loadings(cluster_9_filtered, reduction = "pca"))

# Add gene names as a column
pca_loadings$Gene <- rownames(pca_loadings)

# Filter genes with positive PC1 loadings and sort hierarchically
positive_genes <- pca_loadings %>%
  filter(PC_1 > 0) %>%  # Positive PC1 loadings
  arrange(desc(PC_1), desc(PC_2), desc(PC_3)) %>%  # Hierarchical sorting
  slice_head(n = 100)  # Top 100 genes

# Filter genes with negative PC1 loadings and sort hierarchically
negative_genes <- pca_loadings %>%
  filter(PC_1 < 0) %>%  # Negative PC1 loadings
  arrange(PC_1, PC_2, PC_3) %>%  # Hierarchical sorting (low to high)
  slice_head(n = 100)  # Top 100 genes

# Prepare data for heatmap (positive genes)
positive_heatmap_data <- positive_genes %>%
  select(-Gene) %>%
  as.matrix()
rownames(positive_heatmap_data) <- positive_genes$Gene

# Prepare data for heatmap (negative genes)
negative_heatmap_data <- negative_genes %>%
  select(-Gene) %>%
  as.matrix()
rownames(negative_heatmap_data) <- negative_genes$Gene

# Heatmap for positive PC1 loadings
ComplexHeatmap::Heatmap(
  positive_heatmap_data,
  name = "Loading",
  col = colorRamp2(c(-max(abs(positive_heatmap_data)), 0, max(abs(positive_heatmap_data))), c("blue", "white", "red")),
  row_names_gp = gpar(fontsize = 8),  # Adjust row name font size
  column_names_gp = gpar(fontsize = 12),  # Adjust column name font size
  column_title = "Principal Components (Positive PC1)",
  row_title = "Genes",
  cluster_rows = FALSE,  # Maintain sorting
  cluster_columns = FALSE  # No clustering
)

# Heatmap for negative PC1 loadings
ComplexHeatmap::Heatmap(
  negative_heatmap_data,
  name = "Loading",
  col = colorRamp2(c(-max(abs(negative_heatmap_data)), 0, max(abs(negative_heatmap_data))), c("blue", "white", "red")),
  row_names_gp = gpar(fontsize = 8),  # Adjust row name font size
  column_names_gp = gpar(fontsize = 12),  # Adjust column name font size
  column_title = "Principal Components (Negative PC1)",
  row_title = "Genes",
  cluster_rows = FALSE,  # Maintain sorting
  cluster_columns = FALSE  # No clustering
)


#recalculate loading along principal curve
library(princurve)
library(dplyr)
library(ggplot2)
library(scales)

# Extract PCA coordinates
pca_coords <- Embeddings(cluster_9_filtered, reduction = "pca")[, 1:3]  # PCA 1, 2, and 3

# Fit a principal curve
fit <- principal_curve(as.matrix(pca_coords))

# Extract the arc length (trajectory along the curve)
cluster_9_filtered$arc_length <- fit$lambda

# Regress each gene onto the arc length to compute loadings
gene_expression <- as.matrix(GetAssayData(cluster_9_filtered, assay = "RNA", slot = "data"))
loadings <- apply(gene_expression, 1, function(gene) {
  lm(gene ~ cluster_9_filtered$arc_length)$coefficients[2]  # Slope as the loading
})

# Organize loadings into a data frame
loadings_df <- data.frame(
  Gene = rownames(gene_expression),
  Loading = loadings
)

# Sort genes by loading values
highest_genes <- loadings_df %>% arrange(desc(Loading)) %>% slice_head(n = 100)
lowest_genes <- loadings_df %>% arrange(Loading) %>% slice_head(n = 100)

# Combine top and bottom genes
top_bottom_genes <- bind_rows(
  highest_genes %>% mutate(Category = "Positive"),
  lowest_genes %>% mutate(Category = "Negative")
)

# Dot plot with heatmap color scale
ggplot(top_bottom_genes, aes(
  x = Loading,
  y = reorder(Gene, Loading),
  color = Loading
)) +
  geom_point(size = 3) +  # Fixed size dots
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Loading"
  ) +
  labs(
    title = "Top 100 Positive and Negative Genes by Loading Along Polynomial Fit Axis",
    x = "Loading",
    y = "Gene"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )


#pc variance

library(ggplot2)

# Extract PCA variance explained
pca_variance <- cluster_9_filtered@reductions$pca@stdev^2  # Eigenvalues
total_variance <- sum(pca_variance)  # Total variance
percent_variance <- (pca_variance / total_variance) * 100  # Percentage variance explained

# Prepare data for plotting
variance_df <- data.frame(
  PC = seq_along(percent_variance),
  Variance = percent_variance
)

# Plot percentage variance explained
ggplot(variance_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  labs(
    title = "Percentage Variance Explained by Principal Components",
    x = "Principal Component",
    y = "Variance Explained (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Optionally: Display cumulative variance explained
variance_df$cumulative <- cumsum(variance_df$Variance)

ggplot(variance_df, aes(x = PC, y = cumulative)) +
  geom_line(color = "darkorange", size = 1) +
  geom_point(color = "darkorange", size = 2) +
  labs(
    title = "Cumulative Variance Explained by Principal Components",
    x = "Principal Component",
    y = "Cumulative Variance Explained (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )


#principal curve in N dimensions
library(princurve)
library(ggplot2)

# Assuming `pca_coords` contains N-dimensional PCA coordinates
pca_coords <- Embeddings(cluster_9_filtered, reduction = "pca")[, 1:15]  # Use first 15 PCs

# Fit the principal curve
fit <- principal_curve(as.matrix(pca_coords))

# Extract curve coordinates (projected into 2D)
curve_2d <- as.data.frame(fit$s[, 1:2])
colnames(curve_2d) <- c("PC1", "PC2")

# Add arc length to the data points
pca_2d <- as.data.frame(pca_coords[, 1:2])  # First two PCs for visualization
colnames(pca_2d) <- c("PC1", "PC2")
pca_2d$arc_length <- fit$lambda
# Visualize the data and the principal curve in 2D
# Visualize the data points in 2D
# Minimal script for scatterplot of data points
ggplot() +
  # Plot the data points
  geom_point(
    data = pca_2d,
    aes(x = PC1, y = PC2, color = arc_length),
    size = 1.5,
    alpha = 0.8
  ) +
  # Plot the principal curve as dots
  geom_point(
    data = curve_2d,
    aes(x = PC1, y = PC2),
    color = "black",
    size = 2
  ) +
  # Add color scale and labels
  scale_color_gradient(low = "blue", high = "red", name = "Arc Length") +
  labs(
    title = "Principal Curve in 2D (Dots Only)",
    x = "PC1",
    y = "PC2"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )








#Plot gene by gene in unrolled space

library(ggplot2)
library(dplyr)

# Define output folder for saving plots
output_folder <- "unrolled_gene_scatterplots2"
if (!dir.exists(output_folder)) dir.create(output_folder)

# Extract gene expression matrix
gene_expression <- as.matrix(GetAssayData(cluster_9_filtered, assay = "RNA", slot = "data"))

# Arc length from the principal curve
arc_length <- cluster_9_filtered$arc_length

# Calculate loadings for each gene
principal_curve_loadings <- apply(gene_expression, 1, function(gene) {
  lm(gene ~ arc_length)$coefficients[2]  # Slope as the loading
})

# Convert loadings to a data frame
loadings_df <- data.frame(
  Gene = rownames(gene_expression),
  Loading = principal_curve_loadings
)

# Identify top 20 positive and negative loading genes
top_positive <- loadings_df %>%
  arrange(desc(Loading)) %>%
  slice_head(n = 200)

top_negative <- loadings_df %>%
  arrange(Loading) %>%
  slice_head(n = 200)

# Combine top contributors
top_genes <- bind_rows(
  top_positive %>% mutate(Type = "Positive"),
  top_negative %>% mutate(Type = "Negative")
)

#save top genes to csv
write.csv(top_genes, "top200_pos_neg_loading.csv")

# Subset expression data for the top genes
top_gene_expression <- gene_expression[top_genes$Gene, , drop = FALSE]

# Convert to long format for plotting
expression_long <- as.data.frame(top_gene_expression) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Cell", values_to = "Expression") %>%
  mutate(
    ArcLength = arc_length[match(Cell, colnames(gene_expression))],  # Match arc length
    GeneType = top_genes$Type[match(Gene, top_genes$Gene)]           # Match gene type
  )

# Plot each gene separately and save
unique_genes <- unique(expression_long$Gene)
for (gene in unique_genes) {
  gene_data <- expression_long %>% filter(Gene == gene)
  
  p <- ggplot(gene_data, aes(x = ArcLength, y = Expression, fill = GeneType)) +
    geom_point(shape = 21, size = 0.5, alpha = 0.8, color = "black", stroke = 0.05) +  # Circles with border
    labs(
      title = paste("Expression of", gene, "in Unrolled Space"),
      x = "Arc Length (Unrolled Space)",
      y = "Expression",
      fill = "Gene Type"
    ) +
    scale_fill_manual(values = c("Positive" = "#FE6997", "Negative" = "#4888FD")) +  # Pastel red and blue
    theme_minimal(base_size = 14) +
    theme(
      plot.background = element_rect(fill = "white"),  # Ensure white background
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "right"
    )
  
  # Save the plot
  ggsave(filename = file.path(output_folder, paste0(gene, ".png")), plot = p, width = 6, height = 4)
}

