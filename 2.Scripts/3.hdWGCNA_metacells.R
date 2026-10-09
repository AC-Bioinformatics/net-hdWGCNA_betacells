# hdWGCNA part 1: set up and build metacells on GSE221156 beta cells, ND vs T2D only.
# All cells are beta cells, so there's no cell type to select - the network is built
# on every cell. Metacells average the authors' normalised values (X, the "data"
# layer), so the network uses the authors' SoupX-corrected, log-normalised data.
# Saves the object so 4.hdWGCNA_network.R can build the network after the metacell
# UMAP has been checked.
Root = '/home/alan.culligan/network'
Data_dir = file.path(Root, '1.Data', 'GSE221156')
Out_dir = file.path(Root, '3.Results')
setwd(Root)
dir.create(Out_dir, showWarnings = FALSE, recursive = TRUE)

# single-cell analysis package
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
theme_set(theme_cowplot())

# Random seed for reproducibility
set.seed(24)

# Load the Seurat object
seu = readRDS(file.path(Data_dir, "GSE221156_beta_seurat.rds"))

# Update to v5 object if seurat is v4 or older
seu <- SeuratObject::UpdateSeuratObject(seu)

# PCA on the authors' normalised values comes from 2.loadRDS.R (saved in the .rds) -
# metacells need a reduction to find each cell's nearest neighbours
if (!"pca" %in% Reductions(seu)) {
  stop("No pca reduction - run 2.loadRDS.R first")
}

# Non-diabetic vs T2D only - prediabetic donors are left out
seu <- subset(seu, subset = condition %in% c("ND", "T2D"))
seu$condition <- droplevels(seu$condition)
seu$donor <- droplevels(seu$donor)
print(table(seu$condition))

# Set up the Seurat object for WGCNA analysis
# Requires a fraction of cells that a gene needs to be expressed in order to be included in the analysis.
seu <- SetupForWGCNA(
  seurat_obj = seu,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed to be included
  wgcna_name = "beta" # the name of the hdWGCNA experiment
)
message("hdWGCNA genes: ", length(GetWGCNAGenes(seu)))

# Constructing metacells within each condition (ND, T2D) from the authors' normalised
# values: averaged, not summed, since they are already normalised.
seu <- MetacellsByGroups(
  seurat_obj = seu,
  group.by = "condition", # metacells never mix ND and T2D cells
  reduction = 'pca', # PCA on the authors' normalised values; no batch correction needed
  dims = 1:20, # same 20 PCs as the UMAP in 2.loadRDS.R
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'condition', # set the Idents of the metacell seurat object
  layer = 'data', # authors' normalised values (Seurat v5)
  slot = 'data', # same, for Seurat v4
  mode = 'average' # average, not sum - the values are already normalised
)

# No NormalizeMetacells(): the averaged values are already normalised. hdWGCNA stores
# them in the metacell "counts" layer, so copy them into "data" for SetDatExpr.
metacell_obj <- GetMetacellObject(seu)
LayerData(metacell_obj, layer = "data") <- LayerData(metacell_obj, layer = "counts")
seu <- SetMetacellObject(seu, metacell_obj)
message("Metacells: ", ncol(metacell_obj))
print(table(metacell_obj$condition))

# Processing the metacells: Scale, PCA and UMAP (the metacell object only keeps the
# group.by column, so it's plotted by condition)
seu <- ScaleMetacells(seu, features=GetWGCNAGenes(seu))
seu <- RunPCAMetacells(seu, features=GetWGCNAGenes(seu))
seu <- RunUMAPMetacells(seu, reduction='pca', dims=1:15)

metacell_plot <- DimPlotMetacells(seu, group.by='condition') + umap_theme() + ggtitle("Metacells by condition")
ggsave(file.path(Out_dir, "7.Metacell_UMAP_plot.pdf"), metacell_plot, width = 7, height = 6)

# Save for 4.hdWGCNA_network.R
saveRDS(seu, file.path(Data_dir, "GSE221156_beta_hdwgcna.rds"))
message("Saved: ", file.path(Data_dir, "GSE221156_beta_hdwgcna.rds"))
