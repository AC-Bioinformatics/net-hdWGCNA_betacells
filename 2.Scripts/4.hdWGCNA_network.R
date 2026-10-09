# hdWGCNA part 2: co-expression network on the ND vs T2D beta-cell metacells, using the
# expression matrix and soft-power test from 3.hdWGCNA_metacells.R. Module eigengenes
# are computed on the single cells using the authors' normalised values (the "data"
# layer from 1.read_GSE221156.R).
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

# Soft power for the network: check 8.SoftPower_plot.pdf / 8.SoftPower_table.csv from
# 3.hdWGCNA_metacells.R. NULL = hdWGCNA's choice (lowest power with SFT.R.sq >= 0.8)
soft_power <- 12 # first power with SFT.R.sq >= 0.8 (0.89; mean connectivity 23)

# Load the object from 3.hdWGCNA_metacells.R (metacells, expression matrix, power table)
seu = readRDS(file.path(Data_dir, "GSE221156_beta_hdwgcna.rds"))
print(GetPowerTable(seu))

# construct co-expression network:
seu <- ConstructNetwork(
  seurat_obj = seu,
  soft_power = soft_power, # NULL = chosen automatically from the power table
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
