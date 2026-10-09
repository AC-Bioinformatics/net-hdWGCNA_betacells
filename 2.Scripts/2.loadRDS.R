### Read the seurat object and file structure
Root = ''
Data_dir = file.path(Root, '')
Out_dir = file.path(Root, '')
setwd(Root)

# Libraries
library(Seurat)

# readRDS
seu = readRDS(file.path(Out_dir, 'GSE221156_beta_seurat.rds'))
sink(file.path(Out_dir, "GSE221156_summary.txt"))
colnames(seu@meta.data)
str(seu@meta.data)
sink()



