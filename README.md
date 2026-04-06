# LCNE preprocessing for the snRNAseq data

Flteirng to get your targeted cells. Runt he script based on the order. 

## 1. 01_build_full_object
- load the data and run QC
- run harmony batch correction 
- save as a seaprate rds file `LC_clustered_harmony`
## 2. 02_publication_fig
- create the figs (note this is highly linked to the other 02 script)

## 3. 02_define_LC_populations
- define the LC population from your clusters

## 4. 03_LC_subclustering 
- from the defined LC clusters check if there are outliers and then remove them 

## 5. 04_export_LC_to_h5ad
- save the results as h5ad for the downstream python analysis 