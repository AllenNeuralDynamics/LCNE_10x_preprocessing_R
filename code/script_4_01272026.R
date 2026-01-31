### convert rds into the h5ad files
## reference: `scratch_shuonan/scripts/LC_NE_dataintegration/rds_to_h5.R`



library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)

packageVersion("Seurat") # ‘5.4.0’
packageVersion("SeuratObject") # ‘5.3.0’


happytargetfile = '../results/LC_subclusters_filtered.rds' # output from script 3
outputfile = '../results/snRNAseq_LCNE.h5ad'
seurat_obj <- readRDS(happytargetfile)
seurat_obj <- SeuratObject::UpdateSeuratObject(seurat_obj)


sce <- as.SingleCellExperiment(seurat_obj) # seurat function 
writeH5AD(sce, outputfile) #zellkonverter function 
