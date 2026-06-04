rm(list = ls())
file_loc     <- "/data/LCv3/CS_Neuromodulatory_Noradrenergic_260513/CS_Neuromodulatory_Noradrenergic_260513-counts.h5ad"
# metafile_loc <- "/data/LCv3/CS_Neuromodulatory_Noradrenergic_260513/CS_Neuromodulatory_Noradrenergic_260513.csv"
samp.dat <- read.csv(metafile_loc, stringsAsFactors = FALSE)  # this is not the correct size



library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(anndata)
library(Matrix)
source("/code/func/integration_functions.R")

# --- Load h5ad + CSV, build a Seurat object equivalent to the old mat + samp.dat ---
ad_obj <- read_h5ad(file_loc)

# anndata: cells x genes -> Seurat wants genes x cells
mat <- Matrix::t(ad_obj$X)   # 32285 398912
metadata <- ad_obj$obs
rownames(mat) <- ad_obj$var[['gene_name']]
# 398912 cells to start with 





seurat_raw <- CreateSeuratObject(counts = mat, meta.data = metadata)
rm(ad_obj, mat, samp.dat); gc()



# --- Pipeline ---
mysample <- load_and_create_seurat(seurat_raw, "Sample3",
                                   doublet_scores_col = "doublet_score")
combined <- integrate_samples_hierarchical(mysample)   #29617 231394


saveRDS(combined, "../scratch/LC_clustered_harmony.rds")
print("DONE!")



################ expecter output from above:################
# Cells after study filter: 398912
# Cells after CreateSeuratObject: 398912
# Cells after QC filtering: 231394 (removed 167518)
# Normalizing layer: counts
# Performing log-normalization

