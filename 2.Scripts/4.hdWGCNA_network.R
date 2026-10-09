# hdWGCNA part 2: co-expression network on the ND vs T2D beta-cell metacells built by
# 3.hdWGCNA_metacells.R. Module eigengenes are computed on the single cells using the
# authors' normalised values (the "data" layer from 1.read_GSE221156.R).
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

# Load the object with metacells from 3.hdWGCNA_metacells.R
seu = readRDS(file.path(Data_dir, "GSE221156_beta_hdwgcna.rds"))

# Set expression matrix to use for the analysis: the normalised metacells.
# No group.by / group_name - there's only one cell type, so every metacell (ND and
# T2D) goes into one network and modules can be compared between conditions.
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
ggsave(file.path(Out_dir, "8.SoftPower_plot.png"), soft_power_plot, width = 10, height = 6)

# Table of soft powers and their corresponding scale-free topology fit indices
power_table <- GetPowerTable(seu)
print(power_table)

# construct co-expression network:
seu <- ConstructNetwork(
  seurat_obj = seu,
  tom_name = 'Beta' # name of the topological overlap matrix written to disk
)

png(file.path(Out_dir, "9.Dendrogram_plot.png"), width = 10, height = 6, units = "in", res = 300)
PlotDendrogram(seu, main='Beta hdWGCNA Dendrogram')
dev.off()
TOM <- GetTOM(seu)

modules <- GetModules(seu)
print(table(modules$module))

# ModuleEigengenes(group.by.vars = ...) needs ScaleData to have been run on the object
if (!any(grepl("ScaleData", names(seu@commands)))) {
  seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
}

# compute all MEs in the full single-cell dataset (uses the authors' "data" layer),
# harmonised by chemistry - not donor, since each donor has only one condition and
# correcting donor would also remove the ND vs T2D difference
seu <- ModuleEigengenes(
  seurat_obj = seu,
  group.by.vars = "chemistry"
)

saveRDS(seu, file.path(Data_dir, "GSE221156_beta_hdwgcna.rds"))
message("Saved: ", file.path(Data_dir, "GSE221156_beta_hdwgcna.rds"))
