# hdWGCNA co-expression network on GSE221156 beta cells. All cells are beta cells,
# so there's no cell type to select - the network is built on every cell.
# PCA/Harmony and the module eigengenes use the authors' normalised values (the
# "data" layer from 1.read_GSE221156.R); metacells are built from raw counts and
# re-normalised, as standard.
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
library(harmony)
theme_set(theme_cowplot())

# Random seed for reproducibility
set.seed(24)

# Load the Seurat object
seu = readRDS(file.path(Data_dir, "GSE221156_beta_seurat.rds"))

# Update to v5 object if seurat is v4 or older
seu <- SeuratObject::UpdateSeuratObject(seu)

# PCA + Harmony on the authors' normalised values come from 2.loadRDS.R (saved in the
# .rds) - metacells need a reduction to find each cell's nearest neighbours.
# ScaleData is also required by ModuleEigengenes(group.by.vars = ...) at the end.
if (!"harmony" %in% Reductions(seu)) {
  stop("No harmony reduction - run 2.loadRDS.R first")
}

# Set up the Seurat object for WGCNA analysis
# Requires a fraction of cells that a gene needs to be expressed in order to be included in the analysis.
seu <- SetupForWGCNA(
  seurat_obj = seu,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed to be included
  wgcna_name = "beta" # the name of the hdWGCNA experiment
)

# Constructing metacells within each donor (so no metacell mixes donors).
# All cells are beta cells, so condition + donor replaces the tutorial's cell_type + Sample;
# each donor has one condition, so condition just carries the label onto the metacells.
# Donors with < 100 beta cells are skipped (min_cells).
seu <- MetacellsByGroups(
  seurat_obj = seu,
  group.by = c("condition", "donor"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'harmony', # select the dimensionality reduction to perform KNN on
  dims = 1:20, # same 20 Harmony dimensions as the UMAP in 2.loadRDS.R
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'condition' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seu <- NormalizeMetacells(seu)

# Processing the metacells: Scale and PCA, Harmony integration, and UMAP visualization
seu <- ScaleMetacells(seu, features=GetWGCNAGenes(seu))
seu <- RunPCAMetacells(seu, features=GetWGCNAGenes(seu))
seu <- RunHarmonyMetacells(seu, group.by.vars='donor')
seu <- RunUMAPMetacells(seu, reduction='harmony', dims=1:15)

p1 <- DimPlotMetacells(seu, group.by='condition') + umap_theme() + ggtitle("Condition")
p2 <- DimPlotMetacells(seu, group.by='donor') + umap_theme() + ggtitle("Donor")

metacell_plot = p1 | p2
ggsave(file.path(Out_dir, "6.Metacell_UMAP_plot.png"), metacell_plot, width = 10, height = 5)

# Set expression matrix to use for the analysis: the normalised metacells.
# No group.by / group_name - there's only one cell type, so every metacell (ND, PD
# and T2D) goes into one network and modules can be compared between conditions.
seu <- SetDatExpr(
  seurat_obj = seu,
  assay = 'RNA', # using RNA assay
  layer = 'data' # normalised metacell expression
)

# Test different soft powers:
seu <- TestSoftPowers(
  seurat_obj = seu,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seu)

# assemble with patchwork
soft_power_plot <- wrap_plots(plot_list, ncol=2)
ggsave(file.path(Out_dir, "7.SoftPower_plot.png"), soft_power_plot, width = 10, height = 6)

# Table of soft powers and their corresponding scale-free topology fit indices
power_table <- GetPowerTable(seu)
head(power_table)

# construct co-expression network:
seu <- ConstructNetwork(
  seurat_obj = seu,
  tom_name = 'Beta' # name of the topological overlap matrix written to disk
)

png(file.path(Out_dir, "8.Dendrogram_plot.png"), width = 10, height = 6, units = "in", res = 300)
PlotDendrogram(seu, main='Beta hdWGCNA Dendrogram')
dev.off()
TOM <- GetTOM(seu)

# compute all MEs in the full single-cell dataset (uses the authors' "data" layer),
# harmonised by donor
seu <- ModuleEigengenes(
  seurat_obj = seu,
  group.by.vars = "donor"
)
