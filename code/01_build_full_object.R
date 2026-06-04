install.packages("anndata")
rm(list = ls())
# file_loc <- "../data/LCv2/10xV4_Neuromodulatory_Noradrenergic_complete_RTX-4134_mat_241111.rda"
# metafile_loc <- "../data/LCv2/10xV4_Neuromodulatory_Noradrenergic_complete_RTX-4134_CR8_samp.dat_241111.rda"



library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(anndata)        # for read_h5ad()
source("./func/integration_functions.R")

# --- Load h5ad counts + CSV metadata, then build a Seurat object ---
ad   <- read_h5ad(file_loc)
meta <- read.csv(metafile_loc, row.names = 1, stringsAsFactors = FALSE)

# anndata stores cells as rows, genes as cols; Seurat wants genes x cells
counts <- Matrix::t(ad$X)
rownames(counts) <- ad$var_names
colnames(counts) <- ad$obs_names

# Align metadata rows to the count matrix columns
meta <- meta[colnames(counts), , drop = FALSE]
seurat_obj <- CreateSeuratObject(counts = counts, meta.data = meta)




# Hand off to the existing pipeline
mysample <- load_and_create_seurat(seurat_obj, meta, "Sample3",
                                   doublet_scores_col = "doublet_score")
combined <- integrate_samples_hierarchical(mysample)
saveRDS(combined, "../scratch/LC_clustered_harmony.rds")
print("DONE!")

# 
# mysample <- load_and_create_seurat(file_loc,  metafile_loc, "Sample3",
#                                    doublet_scores_col = "doublet_score")    # 398912 -> 234813 -> 231500 (after QC)
# combined <- integrate_samples_hierarchical(mysample)
# saveRDS(combined, "../scratch/LC_clustered_harmony.rds")
# print("DONE!")
# 
# 
