rm(list = ls())
file_loc <- "../data/LCv2/10xV4_Neuromodulatory_Noradrenergic_complete_RTX-4134_mat_241111.rda"
metafile_loc <- "../data/LCv2/10xV4_Neuromodulatory_Noradrenergic_complete_RTX-4134_CR8_samp.dat_241111.rda"

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
source("./func/integration_functions.R")

mysample <- load_and_create_seurat(file_loc,  metafile_loc, "Sample3",
  doublet_scores_col = "doublet_score")    # 398912 -> 234813 -> 231500 (after QC)
combined <- integrate_samples_hierarchical(mysample)
saveRDS(combined, "../scratch/LC_clustered_harmony_1.rds")
print("DONE!")


