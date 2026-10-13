# hdWGCNA part 5: the co-expression network of one module (brown by default).
#   1. Module genes and their topological overlap (TOM) from 4.hdWGCNA_network.R
#   2. Hub genes and connectivity: kME (module membership, script 5) and TOM-based
#      intramodular connectivity (kWithin), plus degree / betweenness in the thresholded network
#   3. hdWGCNA network plots (ModuleNetworkPlot, HubGeneNetworkPlot, PlotKMEs) and the
#      WGCNA kME vs kWithin scatterplot
#   4. Cytoscape tables: edges (strongest TOM pairs) and node attributes, loaded into
#      Cytoscape by 7b.network_cytoscape.R (RCy3)

# Path configuration
Root <- "/home/alan.culligan/network"
Data_dir <- file.path(Root, "1.Data", "GSE221156")
Out_dir <- file.path(Root, "3.Results")
setwd(Root)

# Libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(WGCNA)
  library(hdWGCNA)
  library(igraph)
})
source(file.path(Root, "2.Scripts", "Functions", "qc_plots.R"))

# HubGeneNetworkPlot uses future.apply, which exports the Seurat object as a global:
# run it sequentially and allow globals the size of the object
future::plan("sequential")
options(future.globals.maxSize = 16 * 1024^3)

# Settings
mod           <- "brown" # module to extract
n_hubs        <- 10      # hub genes = top n_hubs by kME (as in GetHubGenes)
edge_quantile <- 0.95    # Cytoscape edges: TOM pairs above this quantile of within-module TOM

mod_dir <- file.path(Out_dir, "12.Network", mod)
dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
stem <- function(x) file.path(mod_dir, paste0(mod, "_", x))

# Load the object with modules and kME (scripts 4-5)
seu <- readRDS(file.path(Data_dir, "GSE221156_beta_hdwgcna.rds"))
modules <- GetModules(seu)
if (!mod %in% modules$module) stop("No module called '", mod, "'")

## 1. Module genes and TOM
# GetTOM loads the TOM written by ConstructNetwork, with genes in GetModules order
TOM <- GetTOM(seu)
diag(TOM) <- 0
genes <- modules$gene_name[modules$module == mod]
message(mod, " module: ", length(genes), " genes")

## 2. Connectivity and hub genes
# intramodularConnectivity on the TOM: kTotal (all genes), kWithin (own module),
# kOut = kTotal - kWithin, kDiff = kWithin - kOut
conn <- intramodularConnectivity(TOM, as.character(modules$module), scaleByMax = FALSE)
TOM_mod <- TOM[genes, genes]
rm(TOM); invisible(gc())

hub_tab <- modules %>%
  filter(module == mod) %>%
  transmute(gene_name, kME = .data[[paste0("kME_", mod)]]) %>%
  bind_cols(conn[genes, c("kTotal", "kWithin", "kOut", "kDiff")]) %>%
  mutate(kWithin_scaled = kWithin / max(kWithin),
         rank_kME       = rank(-kME, ties.method = "first"),
         rank_kWithin   = rank(-kWithin, ties.method = "first"),
         hub            = if_else(rank_kME <= n_hubs, "hub", "other"))

# Thresholded network: within-module TOM pairs above edge_quantile
tom_upper <- TOM_mod[upper.tri(TOM_mod)]
tom_threshold <- unname(quantile(tom_upper, edge_quantile))
edges <- which(upper.tri(TOM_mod) & TOM_mod > tom_threshold, arr.ind = TRUE)
edges <- data.frame(source = genes[edges[, 1]], target = genes[edges[, 2]],
                    weight = TOM_mod[edges])
g <- graph_from_data_frame(edges, directed = FALSE, vertices = genes)
hub_tab <- hub_tab %>%
  mutate(degree      = degree(g)[gene_name],
         betweenness = betweenness(g, weights = 1 / E(g)$weight)[gene_name])

# Pathway labels from 6.module_pathways.R: the most significant Hallmark term each gene is in
hallmark_csv <- file.path(Out_dir, "11.Module_pathways", paste0(mod, "_ORA_Hallmark"), "tables",
                          paste0(mod, "_ORA_Hallmark_results.csv"))
hub_tab$hallmark <- "none"
if (file.exists(hallmark_csv)) {
  hallmark <- read.csv(hallmark_csv) %>% arrange(p.adjust)
  for (i in rev(seq_len(nrow(hallmark)))) {
    in_term <- hub_tab$gene_name %in% strsplit(hallmark$geneID[i], "/")[[1]]
    hub_tab$hallmark[in_term] <- sub("^HALLMARK_", "", hallmark$ID[i])
  }
}

hub_tab <- arrange(hub_tab, rank_kME)
write.csv(hub_tab, stem("node_attributes.csv"), row.names = FALSE)
print(head(hub_tab, n_hubs))

## 3. Plots
# hdWGCNA: network of the top hub genes (writes <outdir>/<module>.pdf)
ModuleNetworkPlot(seu, mods = mod, outdir = mod_dir, n_inner = n_hubs)

pdf(stem("HubGeneNetworkPlot.pdf"), width = 10, height = 10)
hub_plot <- HubGeneNetworkPlot(seu, mods = mod, n_hubs = n_hubs, n_other = 20, edge_prop = 0.75)
if (inherits(hub_plot, "gg")) print(hub_plot)
dev.off()

pdf(stem("PlotKMEs.pdf"), width = 14, height = 10)
print(PlotKMEs(seu, n_hubs = n_hubs, ncol = 4))
dev.off()

# WGCNA: module membership vs intramodular connectivity
pdf(stem("kME_vs_kWithin.pdf"), width = 7, height = 7)
verboseScatterplot(hub_tab$kME, hub_tab$kWithin_scaled, col = mod,
                   main = paste0(mod, " module: kME vs scaled kWithin\n"),
                   xlab = paste0("kME (", mod, ")"), ylab = "kWithin / max(kWithin)")
dev.off()

# WGCNA scale-free check on intramodular connectivity (kWithin), plus the degree
# distribution of the thresholded network
sft <- scaleFreeFitIndex(hub_tab$kWithin, nBreaks = 10)
pdf(stem("scale_free.pdf"), width = 15, height = 5)
par(mfrow = c(1, 3))
hist(hub_tab$kWithin, breaks = 30, col = mod,
     main = paste0(mod, " module: kWithin distribution"), xlab = "kWithin (TOM)")
scaleFreePlot(hub_tab$kWithin, nBreaks = 10, main = paste0(mod, " module: scale-free fit (kWithin)\n"))
hist(hub_tab$degree, breaks = 30, col = mod,
     main = paste0(mod, " module: degree (TOM > ", signif(tom_threshold, 3), ")"),
     xlab = "Degree")
dev.off()

## 4. Cytoscape tables
# Edge and node files in WGCNA's Cytoscape format (nodes = genes with >= 1 edge)
exportNetworkToCytoscape(TOM_mod,
                         edgeFile  = stem("cytoscape_edges.txt"),
                         nodeFile  = stem("cytoscape_nodes.txt"),
                         weighted  = TRUE,
                         threshold = tom_threshold,
                         nodeNames = genes,
                         nodeAttr  = rep(mod, length(genes)))

# Summary txt
write_network_summary(
  stem("network_summary.txt"), hub_tab,
  params = list(
    "Module"             = mod,
    "Object"             = file.path(Data_dir, "GSE221156_beta_hdwgcna.rds"),
    "Hub genes"          = paste("top", n_hubs, "by kME"),
    "Edge quantile"      = edge_quantile,
    "TOM threshold"      = signif(tom_threshold, 4),
    "Scale-free R^2"     = signif(sft$Rsquared.SFT, 3),
    "Scale-free slope"   = signif(sft$slope.SFT, 3),
    "Trunc. exp. adj R^2" = signif(sft$truncatedExponentialAdjRsquared, 3)
  ),
  n_edges  = nrow(edges),
  n_pairs  = length(tom_upper),
  files    = list("Node attributes" = stem("node_attributes.csv"),
                  "Cytoscape edges" = stem("cytoscape_edges.txt"),
                  "Cytoscape nodes" = stem("cytoscape_nodes.txt"),
                  "Plots"           = mod_dir)
)
message("Saved: ", mod_dir)
