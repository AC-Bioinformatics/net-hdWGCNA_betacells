# hdWGCNA part 1: set up and build metacells on GSE221156 beta cells, ND vs T2D only.
# All cells are beta cells, so there's no cell type to select - the network is built
# on every cell. Metacells are built from raw counts and re-normalised, as standard.
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
if (!"harmony" %in% Reductions(seu)) {
  stop("No harmony reduction - run 2.loadRDS.R first")
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

# Constructing metacells within each donor (so no metacell mixes donors).
# All cells are beta cells, so condition + donor replaces the tutorial's cell_type + Sample;
# each donor has one condition, so condition just carries the label onto the metacells.
# Donors with < 100 beta cells are skipped (min_cells).
seu <- MetacellsByGroups(
  seurat_obj = seu,
  group.by = c("condition", "donor", "chemistry"), # chemistry is constant within a donor; carried for the plots
  reduction = 'harmony', # select the dimensionality reduction to perform KNN on
  dims = 1:20, # same 20 Harmony dimensions as the UMAP in 2.loadRDS.R
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'condition' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seu <- NormalizeMetacells(seu)
metacell_obj <- GetMetacellObject(seu)
message("Metacells: ", ncol(metacell_obj))
print(table(metacell_obj$condition))

# Processing the metacells: Scale and PCA, then UMAP before and after Harmony.
# Harmony corrects chemistry, not donor - each donor has only one condition, so
# correcting donor would also remove the ND vs T2D difference.
seu <- ScaleMetacells(seu, features=GetWGCNAGenes(seu))
seu <- RunPCAMetacells(seu, features=GetWGCNAGenes(seu))

seu <- RunUMAPMetacells(seu, reduction='pca', dims=1:15)
p1 <- DimPlotMetacells(seu, group.by='condition') + umap_theme() + ggtitle("Condition (uncorrected)")
p2 <- DimPlotMetacells(seu, group.by='chemistry') + umap_theme() + ggtitle("Chemistry (uncorrected)")
p3 <- DimPlotMetacells(seu, group.by='donor') + umap_theme() + NoLegend() + ggtitle("Donor (uncorrected)")

seu <- RunHarmonyMetacells(seu, group.by.vars='chemistry')
seu <- RunUMAPMetacells(seu, reduction='harmony', dims=1:15)
p4 <- DimPlotMetacells(seu, group.by='condition') + umap_theme() + ggtitle("Condition (Harmony: chemistry)")
p5 <- DimPlotMetacells(seu, group.by='chemistry') + umap_theme() + ggtitle("Chemistry (Harmony: chemistry)")
p6 <- DimPlotMetacells(seu, group.by='donor') + umap_theme() + NoLegend() + ggtitle("Donor (Harmony: chemistry)")

metacell_plot <- (p1 | p2 | p3) / (p4 | p5 | p6)
ggsave(file.path(Out_dir, "7.Metacell_UMAP_plot.png"), metacell_plot, width = 18, height = 10)

# Save for 4.hdWGCNA_network.R
saveRDS(seu, file.path(Data_dir, "GSE221156_beta_hdwgcna.rds"))
message("Saved: ", file.path(Data_dir, "GSE221156_beta_hdwgcna.rds"))
