# hdWGCNA part 1: set up, build metacells and test soft powers on GSE221156 beta
# cells, ND vs T2D only, v3 chemistry only.
# All cells are beta cells, so there's no cell type to select - the network is built
# on every cell. Metacells average the authors' normalised values (X, the "data"
# layer), so the network uses the authors' SoupX-corrected, log-normalised data.
# Saves the object so 4.hdWGCNA_network.R can build the network after the metacell
# UMAP and soft-power plots have been checked.
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

# Harmony embedding (as in the authors' methods) comes from 2.loadRDS.R (saved in the
# .rds) - metacells need a reduction to find each cell's nearest neighbours
if (!"harmony" %in% Reductions(seu)) {
  stop("No harmony reduction - run 2.loadRDS.R first")
}

# Non-diabetic vs T2D only - prediabetic donors are left out.
# v3 chemistry only: v2 and v3 separate clearly (UMAP, metacells) and dominated the
# soft-power test (positive slope, mean connectivity ~6,000 at power 1). v2 is also
# unevenly split (31% of T2D vs 17% of ND cells), so chemistry and disease can't be
# separated in a mixed network. The v2 cells can be used later to check whether the
# v3 modules change in the same direction between ND and T2D.
seu <- subset(seu, subset = condition %in% c("ND", "T2D") & chemistry == "v3")
seu$condition <- droplevels(seu$condition)
seu$donor <- droplevels(seu$donor)
print(table(seu$condition))
message("Donors per condition (v3): ",
        paste(names(table(unique(seu[[c("donor", "condition")]])$condition)),
              table(unique(seu[[c("donor", "condition")]])$condition), collapse = ", "))

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
  reduction = 'harmony', # Harmony embedding, as in the authors' methods
  dims = 1:20, # same 20 Harmony dimensions as the UMAP (authors' and ours)
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

# Set expression matrix to use for the analysis: the metacells' averaged authors' normalised values (X).
# No group.by / group_name - there's only one cell type, so every metacell (ND and
# T2D) goes into one network and modules can be compared between conditions.
seu <- SetDatExpr(
  seurat_obj = seu,
  assay = 'RNA', # using RNA assay
  layer = 'data' # authors' normalised values, averaged per metacell
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
ggsave(file.path(Out_dir, "8.SoftPower_plot.pdf"), soft_power_plot, width = 10, height = 6)

# Table of soft powers and their corresponding scale-free topology fit indices
power_table <- GetPowerTable(seu)
print(power_table)
write.csv(power_table, file.path(Out_dir, "8.SoftPower_table.csv"), row.names = FALSE)
message("Lowest power with SFT.R.sq >= 0.8: ",
        min(power_table$Power[power_table$SFT.R.sq >= 0.8], na.rm = TRUE))

# Save for 4.hdWGCNA_network.R
saveRDS(seu, file.path(Data_dir, "GSE221156_beta_hdwgcna.rds"))
message("Saved: ", file.path(Data_dir, "GSE221156_beta_hdwgcna.rds"))
