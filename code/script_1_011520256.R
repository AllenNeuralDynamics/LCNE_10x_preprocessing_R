
packageVersion("Seurat") # ‘5.4.0’
packageVersion("SeuratObject") # ‘5.3.0’



if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager", ask = FALSE)

BiocManager::install("SingleCellExperiment", ask = FALSE, update = FALSE)
library(SingleCellExperiment)

if (!requireNamespace("zellkonverter", quietly = TRUE)) {
  BiocManager::install("zellkonverter", ask = FALSE, update = FALSE)
}
library(zellkonverter) # 1.16




# cleaned up version 
rm(list = ls())
gc()  # Garbage collection at start

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------------------------
# Memory-optimized Seurat object creation
# ------------------------------------------------------------------------------

load_and_create_seurat <- function(data_file, metadata_file, sample_name, 
                                   min_genes = 2000, max_genes = 15000, 
                                   max_mt_pct = 4, max_rb_pct = 4, 
                                   doublet_scores_col = 'doublet_score') {
  
  temp_env <- new.env()
  load(data_file, envir = temp_env)
  load(metadata_file, envir = temp_env)
  
  mat <- temp_env$mat
  samp.dat <- temp_env$samp.dat
  
  # Filter metadata
  if (!"studies" %in% colnames(samp.dat)) {
    stop("The metadata file does not contain a 'studies' column.")
  }
  
  samp.dat <- samp.dat[samp.dat$studies == 'Neuromodulatory_Noradrenergic', ]
  
  # Filter matrix and ensure alignment
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
  
  # Clear large objects immediately
  rm(mat, samp.dat, temp_env)
  gc()
  
  # Add metadata
  seurat_object <- AddMetaData(seurat_object, metadata = sample_name, col.name = 'experiment')
  seurat_object <- AddMetaData(seurat_object, metadata = '10xV4', col.name = 'platform')
  
  # Calculate QC metrics
  seurat_object$percent.mt <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  seurat_object$percent.rb <- PercentageFeatureSet(seurat_object, pattern = "^RP[SL]")
  
  # Add doublet scores
  if (!is.null(doublet_scores_col) && doublet_scores_col %in% colnames(seurat_object@meta.data)) {
    # Already added via meta.data in CreateSeuratObject
  } else {
    warning(paste0("Doublet score column '", doublet_scores_col, "' not found in metadata."))
  }
  
  # Filter cells
  if (!is.null(doublet_scores_col) && "doublet_score" %in% colnames(seurat_object@meta.data)) {
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
  }
  
  gc()  # Clean up after subsetting
  
  # Normalize data
  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Find variable features
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 1500, assay = "RNA")
  
  return(seurat_object)
}

# ------------------------------------------------------------------------------
# Memory-optimized Hierarchical Harmony integration
# ------------------------------------------------------------------------------

integrate_samples_hierarchical <- function(seurat_object, 
                                           file_col = "batch_vendor_name", 
                                           port_well_col = "rna_amplification") {
  ## ---- Metadata cleanup ----
  md <- seurat_object@meta.data
  for (col in colnames(md)) {
    md[[col]][is.na(md[[col]])] <- "not_available"
  }
  seurat_object@meta.data <- md
  rm(md); gc()
  
  ## ---- PCA ----
  seurat_object <- ScaleData(
    seurat_object,
    features = VariableFeatures(seurat_object),
    assay = "RNA"
  )
  
  seurat_object <- RunPCA(
    seurat_object,
    features = VariableFeatures(seurat_object),
    assay = "RNA",
    npcs = 30
  )
  
  ## ---- Harmony level 1: vendor / RTX ----
  seurat_object <- RunHarmony(
    object = seurat_object,
    group.by.vars = file_col,
    reduction.use = "pca",
    dims.use = 1:30,
    assay.use = "RNA",
    project.dim = FALSE
  )
  seurat_object@reductions$harmony_batch1 <- seurat_object@reductions$harmony
  
  ## ---- Harmony level 2: batch ----
  if (!"batch" %in% colnames(seurat_object@meta.data)) {
    stop("The metadata column 'batch' does not exist.")
  }
  
  seurat_object <- RunHarmony(
    object = seurat_object,
    group.by.vars = "batch",
    reduction.use = "harmony_batch1",
    dims.use = 1:30
  )
  seurat_object@reductions$harmony_batch2 <- seurat_object@reductions$harmony
  
  ## ---- Harmony level 3: RNA amplification ----
  seurat_object <- RunHarmony(
    object = seurat_object,
    group.by.vars = port_well_col,
    reduction.use = "harmony_batch2",
    dims.use = 1:30
  )
  seurat_object@reductions$harmony_final <- seurat_object@reductions$harmony
  gc()
  
  ## ---- UMAP / graph / clustering ----
  seurat_object <- RunUMAP(
    seurat_object,
    reduction = "harmony_final",
    dims = 1:30
  )
  
  seurat_object <- FindNeighbors(
    seurat_object,
    reduction = "harmony_final",
    dims = 1:30
  )
  
  seurat_object <- FindClusters(
    seurat_object,
    resolution = 0.5,
    random.seed = 42
  )
  
  return(seurat_object)
}

# ------------------------------------------------------------------------------
# Main execution with memory monitoring
# ------------------------------------------------------------------------------
# Monitor memory before starting
print(paste("Memory before loading:", format(object.size(ls()), units = "MB")))
file3 <- "../data/LCv2/10xV4_Neuromodulatory_Noradrenergic_complete_RTX-4134_mat_241111.rda"
metafile3 <- "../data/LCv2/10xV4_Neuromodulatory_Noradrenergic_complete_RTX-4134_CR8_samp.dat_241111.rda"
sample3 <- load_and_create_seurat(
  file3,
  metafile3,
  "Sample3",
  doublet_scores_col = "doublet_score"
)

print(sample3)


print(paste("Memory after loading:", format(object.size(sample3), units = "MB")))
combined <- integrate_samples_hierarchical(sample3)
dim(combined)


rm(sample3)
gc()

### gc()
# used    (Mb)  gc trigger    (Mb)    max used     (Mb)
# Ncells    4514500   241.2     9427370   503.5     9427370    503.5
# Vcells 3803729741 29020.2 11355150692 86633.0 14118705710 107717.2
###


print(paste("Memory before saving:", format(object.size(combined), units = "MB")))
saveRDS(combined, "../scratch/LC_clustered_harmony.rds")
print("Analysis complete!")


