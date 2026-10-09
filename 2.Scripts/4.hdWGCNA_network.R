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
theme_set(theme_cowplot())

# Random seed for reproducibility
set.seed(24)

# Load the object with metacells from 3.hdWGCNA_metacells.R
seu = readRDS(file.path(Data_dir, "GSE221156_beta_hdwgcna.rds"))

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

# construct co-expression network:
seu <- ConstructNetwork(
  seurat_obj = seu,
  tom_name = 'Beta' # name of the topological overlap matrix written to disk
)

pdf(file.path(Out_dir, "9.Dendrogram_plot.pdf"), width = 10, height = 6)
PlotDendrogram(seu, main='Beta hdWGCNA Dendrogram')
dev.off()
TOM <- GetTOM(seu)

modules <- GetModules(seu)
print(table(modules$module))

# compute all MEs in the full single-cell dataset (uses the authors' "data" layer).
# No batch correction here: chemistry goes in as a covariate when MEs are compared
# between ND and T2D (per donor), e.g. ME ~ condition + chemistry
seu <- ModuleEigengenes(
  seurat_obj = seu
)

saveRDS(seu, file.path(Data_dir, "GSE221156_beta_hdwgcna.rds"))
message("Saved: ", file.path(Data_dir, "GSE221156_beta_hdwgcna.rds"))
