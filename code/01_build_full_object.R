rm(list = ls())
file_loc <- "../data/LCv2/10xV4_Neuromodulatory_Noradrenergic_complete_RTX-4134_mat_241111.rda"
metafile_loc <- "../data/LCv2/10xV4_Neuromodulatory_Noradrenergic_complete_RTX-4134_CR8_samp.dat_241111.rda"


library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
source("./func/integration_functions.R")

mysample <- load_and_create_seurat(file_loc,  metafile_loc, "Sample3", doublet_scores_col = "doublet_score")    
# Cells after study filter: 398912
# Cells after CreateSeuratObject: 398912
# Cells after QC filtering: 231500 (removed 167412)
# Cells after normalization: 231500
# Cells after HVG selection: 231500
combined <- integrate_samples_hierarchical(mysample)

### Number of communities: 63
## 7 singletons identified. 56 final clusters.

saveRDS(combined, "../scratch/LC_clustered_harmony.rds")
print("DONE!")  




