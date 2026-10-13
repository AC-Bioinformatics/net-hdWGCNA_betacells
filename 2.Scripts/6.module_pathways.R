# Pathway analysis on modules of interest that were identified in the previous
# script '5.hdWGCNA_module_correlations.R'. Pathway analysis will be conducted
# using ClusterProfiler

# Path configuration
Root <- "/home/alan.culligan/network"
Data_dir <- file.path(Root, "1.Data", "GSE221156")
Out_dir <- file.path(Root, "3.Results")
setwd(Root)

# Libraries
library(clusterProfiler)
library(Seurat)
library(hdWGCNA)

# Load data (seurat object with hdWGCNA modules and eigengenes from scripts 4-5)
seu <- readRDS(file.path(Data_dir, "GSE221156_beta_hdwgcna.rds"))
