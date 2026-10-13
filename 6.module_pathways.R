# Pathway analysis on modules of interest that were identified in the previous
# script '5.hdWGCNA_module_correlations.R. Pathway analysis will be conducted
# using ClusterProfiler

# Path configuration
Root <-
  # Libraries
  library(clusterProfiler)
library(Seurat)

# Load data
seu <- readRDS(file.path(Data_dir, ''))
