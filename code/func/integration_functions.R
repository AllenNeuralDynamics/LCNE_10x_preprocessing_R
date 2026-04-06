
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------------------------
# load and fiter 
# ------------------------------------------------------------------------------


load_and_create_seurat <- function(data_file, metadata_file, sample_name, 
                                   min_genes = 2000, max_genes = 15000, 
                                   max_mt_pct = 4, max_rb_pct = 4, 
                                   doublet_scores_col = "doublet_score",
                                   verbose = TRUE) {
  
  temp_env <- new.env()
  load(data_file, envir = temp_env)
  load(metadata_file, envir = temp_env)
  
  mat <- temp_env$mat
  samp.dat <- temp_env$samp.dat
  
  if (!"studies" %in% colnames(samp.dat)) {
    stop("Missing 'studies' column in metadata.")
  }
  
  samp.dat <- samp.dat[samp.dat$studies == "Neuromodulatory_Noradrenergic", ]
  mat <- mat[, colnames(mat) %in% samp.dat$sample_id]
  
  if (ncol(mat) != nrow(samp.dat)) {
    stop("Counts matrix and metadata mismatch.")
  }
  
  if (verbose) {
    message("Cells after study filter: ", ncol(mat))
  }
  
  seurat_object <- CreateSeuratObject(
    counts = mat,
    meta.data = samp.dat,
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )
  
  if (verbose) {
    message("Cells after CreateSeuratObject: ", ncol(seurat_object))
  }
  
  rm(mat, samp.dat, temp_env); gc()
  
  seurat_object <- AddMetaData(seurat_object, metadata = sample_name, col.name = "experiment")
  seurat_object <- AddMetaData(seurat_object, metadata = "10xV4", col.name = "platform")
  
  seurat_object$percent.mt <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
  seurat_object$percent.rb <- PercentageFeatureSet(seurat_object, pattern = "^rp[sl]")
  seurat_object$percent.mt <- as.numeric(seurat_object$percent.mt)
  hist(seurat_object$percent.mt, breaks = 50, xlab = "percent.mt", col = "gray")
  if (!is.null(doublet_scores_col) &&
      !doublet_scores_col %in% colnames(seurat_object@meta.data)) {
    warning(paste0("Doublet score column '", doublet_scores_col, "' not found."))
  }
  
  n_before <- ncol(seurat_object)
  
  if (!is.null(doublet_scores_col) &&
      doublet_scores_col %in% colnames(seurat_object@meta.data)) {
    
    seurat_object <- subset(
      seurat_object,
      subset =
        nFeature_RNA > min_genes &
        nFeature_RNA < max_genes &
        percent.mt < max_mt_pct &
        percent.rb < max_rb_pct &
        doublet_score < 0.4
    )
  } else {
    seurat_object <- subset(
      seurat_object,
      subset =
        nFeature_RNA > min_genes &
        nFeature_RNA < max_genes &
        percent.mt < max_mt_pct &
        percent.rb < max_rb_pct
    )
  }
  
  if (verbose) {
    message(
      "Cells after QC filtering: ",
      ncol(seurat_object),
      " (removed ", n_before - ncol(seurat_object), ")"
    )
  }
  
  gc()
  
  seurat_object <- NormalizeData(
    seurat_object,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
  
  if (verbose) {
    message("Cells after normalization: ", ncol(seurat_object))
  }
  
  seurat_object <- FindVariableFeatures(
    seurat_object,
    selection.method = "vst",
    nfeatures = 1500,
    assay = "RNA"
  )
  
  if (verbose) {
    message("Cells after HVG selection: ", ncol(seurat_object))
  }
  
  seurat_object
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
# Plotting
# ------------------------------------------------------------------------------

plot_integration_results <- function(combined) {
  p1 <- DimPlot(combined, reduction = "pca", group.by = "orig.ident", pt.size = 0.1) +
    ggtitle("Before Harmony Integration")
  p2 <- DimPlot(combined, reduction = "harmony", group.by = "orig.ident", pt.size = 0.1) +
    ggtitle("After Harmony Integration")
  p3 <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 0.1) +
    ggtitle("Clusters")
  return(list(batch_effects = p1 + p2, clusters = p3))
}
