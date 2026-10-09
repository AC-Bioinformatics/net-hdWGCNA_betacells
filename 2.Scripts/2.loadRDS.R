### Read the seurat object and file structure
Root = '/home/alan.culligan/network'
Data_dir = file.path(Root, '1.Data', 'GSE221156')
Out_dir = file.path(Root, '3.Results')
setwd(Root)

# Libraries
library(Seurat)
library(ggplot2)
source(file.path(Root, "2.Scripts", "Functions", "qc_plots.R"))

# readRDS and look at file structure
summary_file <- file.path(Out_dir, "GSE221156_summary.txt")

# Add a line of text to the end of the file (creates the file if it doesn't exist)
write_log <- function(...) {
  cat(paste(...), "\n", file = summary_file, append = TRUE, sep = "")
}

# Add an object's printed output (data frames, str(), Seurat objects, etc.)
write_obj <- function(x) {
  cat(capture.output(x), file = summary_file, append = TRUE, sep = "\n")
}
seu <- readRDS(file.path(Data_dir, "GSE221156_beta_seurat.rds"))

write_log("\n=====", format(Sys.time()), "=====")
write_log("Seurat object summary:")
write_obj(seu)
write_log("meta.data columns:")
write_log(paste(colnames(seu@meta.data), collapse = ", "))
write_log("meta.data structure:")
write_obj(str(seu@meta.data))
write_log("Assays:", paste(Assays(seu), collapse = ", "))
write_log("DefaultAssay:", DefaultAssay(seu))

# Quick qc check of the Seurat object
write_log("\nCells per condition:")
write_obj(table(seu$condition))
write_log("\nDonors per condition:")
write_obj(tapply(seu$donor, seu$condition, function(x) length(unique(x))))
write_log("\nBeta cells per donor:")
write_obj(table(seu$donor, seu$condition)[, c("ND", "PD", "T2D")])
write_log("\nChemistry by condition:")
write_obj(table(seu$chemistry, seu$condition))

# Mitochondrial, ribosomal genes and hemoglobin genes
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RPL|^RPS")   
seu[["percent.hb"]] <- PercentageFeatureSet(seu, pattern = "^HBA|^HBB|^HBZ")

# Violin plot of QC metrics for each condition (ND, PD, T2D)
vln_plot <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA",
                                      "percent.mt", "percent.ribo", "percent.hb"),
                    group.by = "condition", pt.size = 0.1, ncol = 3, raster = TRUE)

ggsave(file.path(Out_dir, "1.QC_violin_plot.pdf"), vln_plot, width = 10, height = 6)

# QC metric distributions per condition - no filtering, just to look at them.
# To filter later: mad_filters = list(nFeature_RNA = 5, percent.mt = 5) if needed
qc_plots(seu,
         group_by = "condition",
         mad_filters = list(),
         out_dir = Out_dir,
         pdf_name = "2.QC_histograms.pdf",
         log_fun = write_log)

# Our own PCA / Harmony / UMAP on the authors' normalised values (data layer), following
# their beta-cell reintegration: their variable genes -> 100 PCs -> Harmony on sex,
# chemistry and ancestry -> UMAP from the first 20 Harmony dimensions.
# Their UMAP is kept alongside as "umap_authors" for comparison.
if (!"pca" %in% Reductions(seu)) {
  # authors' variable genes are set in 1.read_GSE221156.R; only recalculate if missing
  if (length(VariableFeatures(seu)) == 0) {
    seu <- FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE) # uses vst on raw counts by default
  }
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = 100, verbose = FALSE)
}

if (!"harmony" %in% Reductions(seu)) {
  seu <- harmony::RunHarmony(seu, group.by.vars = c("sex", "chemistry", "self_reported_ethnicity"),
                             verbose = FALSE)
}

# plot pca and explained variance
pca_plot <- DimPlot(seu, reduction = "pca", group.by = "donor", raster = TRUE) +
  ggtitle("PCA by donor")
ggsave(file.path(Out_dir, "3.PCA_plot.pdf"), pca_plot, width = 6, height = 4)

# Elbow plot to determine the number of PCs to use for downstream analysis
elbow_plot <- ElbowPlot(seu, ndims = 100)
ggsave(file.path(Out_dir, "4.Elbow_plot.pdf"), elbow_plot, width = 6, height = 4)

if (!"umap" %in% Reductions(seu)) {
  seu <- RunUMAP(seu, reduction = "harmony", dims = 1:20, verbose = FALSE) # 20 dims, as the authors
}

# plot umap coloured by condition (ND vs T2D) - ours next to the authors'
nd_t2d <- colnames(seu)[seu$condition %in% c("ND", "T2D")]
umap_plot <- DimPlot(seu, reduction = "umap", group.by = "condition", raster = TRUE, cells = nd_t2d) +
  ggtitle("Our UMAP by condition") |
  DimPlot(seu, reduction = "umap_authors", group.by = "condition", raster = TRUE, cells = nd_t2d) +
  ggtitle("Authors' UMAP by condition")
ggsave(file.path(Out_dir, "5.UMAP_plot.pdf"), umap_plot, width = 12, height = 5)

# batch check - after Harmony, cells shouldn't separate by chemistry or donor
batch_plot <- DimPlot(seu, reduction = "umap", group.by = "chemistry", raster = TRUE) +
  DimPlot(seu, reduction = "umap", group.by = "donor", raster = TRUE) + NoLegend()
ggsave(file.path(Out_dir, "6.UMAP_batch_check.pdf"), batch_plot, width = 12, height = 5)

# Save the Seurat object with QC metrics and PCA/UMAP embeddings to override the previous version
saveRDS(seu, file.path(Data_dir, "GSE221156_beta_seurat.rds"))









