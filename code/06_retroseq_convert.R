### convert from RData to something python readable
rm(list = ls())
filepath <- '/data/LCNE_smartseq_raw_NeMO/NeMO_LC_260514//'
filelist = list.files(filepath, pattern="")
metafile = filelist[grep('samp.dat',filelist)]
load(paste0(filepath, metafile))  # samp.dat2
exonfile = filelist[grep('exon',filelist)]  # exon2
intronfile = filelist[grep('intron',filelist)]  # intron2
## load and confirm the dimensions
load(paste0(filepath, intronfile))
dim(intron2)   # 1383 
load(paste0(filepath, exonfile))
dim(exon2) # 32245  1383
allcounts = t(exon2)+t(intron2)
max(allcounts)  # 59784


#### save the file
saveloc = '../results/processed_data/retroseqdata_raw_from_R/' # output from script 3
if (!dir.exists(saveloc)) {
  dir.create(saveloc, recursive = TRUE)
}


# metadata
savedfilename = paste0(saveloc, 'metadata.csv')
write.csv(samp.dat2, savedfilename, row.names = TRUE)

## all counrs (intron+exons)
mat <- as(as.matrix(allcounts), "dgCMatrix")
savedfilename = paste0(saveloc, 'retro_raw.mtx')
writeMM(mat, savedfilename)
write.table(colnames(allcounts), file = paste0(saveloc, 'genes.tsv'),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(rownames(allcounts), file = paste0(saveloc, 'cells.tsv'),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

