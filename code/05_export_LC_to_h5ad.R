### convert rds into the h5ad files
## reference: `scratch_shuonan/scripts/LC_NE_dataintegration/rds_to_h5.R`

library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)
source("./func/integration_functions.R")


packageVersion("Seurat") # ‘5.4.0’
packageVersion("SeuratObject") # ‘5.3.0’


happytargetfile = '../scratch/LC_subclusters_filtered.rds' # output from script 3
outputfile = '../results/processed_data/snRNAseq_LCNE.h5ad'
seurat_obj <- readRDS(happytargetfile)
seurat_obj <- SeuratObject::UpdateSeuratObject(seurat_obj)


sce <- as.SingleCellExperiment(seurat_obj) # seurat function 
dir.create('../results/processed_data/', showWarnings = FALSE)
writeH5AD(sce, outputfile) #zellkonverter function 
