### Script to read the GSE221156 beta-cell h5ad (CELLxGENE, Bandesh et al. 2026)
### into a Seurat object. Uses the raw UMI counts in raw.X (not the authors'
### normalised X) and keeps all cell metadata plus the authors' UMAP.
library(Seurat)
library(rhdf5)
library(Matrix)

# Path configuration
Root <- "/home/alan.culligan/network/"
setwd(Root)

h5ad_file <- "1.Data/GSE221156/GSE221156_beta_cellxgene.h5ad"
rds_out   <- "1.Data/GSE221156/GSE221156_beta_seurat.rds"

h5_names <- h5ls(h5ad_file, recursive = TRUE)
h5_names <- sub("^/", "", file.path(h5_names$group, h5_names$name))

## Read one obs/var column; categoricals are stored as 0-based codes + categories
read_column <- function(group, col) {
  path <- paste0(group, "/", col)
  if (paste0(path, "/codes") %in% h5_names) {
    codes <- as.integer(h5read(h5ad_file, paste0(path, "/codes")))
    cats  <- as.character(h5read(h5ad_file, paste0(path, "/categories")))
    codes[codes < 0] <- NA
    factor(cats[codes + 1], levels = cats)
  } else {
    as.vector(h5read(h5ad_file, path))
  }
}

## 1. Raw counts: CSR (cells x genes) is the same layout as CSC (genes x cells)
shape <- h5readAttributes(h5ad_file, "raw/X")$shape
n_cells <- shape[1]
n_genes <- shape[2]
counts <- new("dgCMatrix",
              i   = as.integer(h5read(h5ad_file, "raw/X/indices")),
              p   = as.integer(h5read(h5ad_file, "raw/X/indptr")),
              x   = as.numeric(h5read(h5ad_file, "raw/X/data")),
              Dim = c(as.integer(n_genes), as.integer(n_cells)))
message("Raw counts: ", n_genes, " genes x ", n_cells, " cells")
stopifnot(all(counts@x == round(counts@x)))

## 2. Gene names: Ensembl IDs in the index, symbols in var$feature_name
raw_ids <- as.character(h5read(h5ad_file, "raw/var/_index"))
var_ids <- as.character(h5read(h5ad_file, "var/_index"))
symbols <- as.character(read_column("var", "feature_name"))[match(raw_ids, var_ids)]
n_unmapped <- sum(is.na(symbols))
gene_names <- ifelse(is.na(symbols), raw_ids, symbols)
gene_names <- gsub("_", "-", gene_names)   # Seurat does not allow underscores
n_dup_symbols <- sum(duplicated(gene_names))
gene_names <- make.unique(gene_names)
message("Genes without a symbol: ", n_unmapped,
        "; duplicated symbols made unique: ", n_dup_symbols)

cell_ids <- as.character(h5read(h5ad_file, "obs/_index"))
dimnames(counts) <- list(gene_names, cell_ids)

## 3. Cell metadata (drop ontology IDs, constant columns and the authors'
##    nCount/nFeature, which Seurat recomputes from the same counts)
obs_cols <- setdiff(h5readAttributes(h5ad_file, "obs")[["column-order"]],
                    c("nCount_RNA", "nFeature_RNA", "observation_joinid",
                      "is_primary_data", "suspension_type", "tissue",
                      "tissue_type", "cell_type"))
obs_cols <- obs_cols[!grepl("_ontology_term_id$", obs_cols)]
meta <- as.data.frame(lapply(setNames(obs_cols, obs_cols),
                             function(col) read_column("obs", col)),
                      check.names = FALSE)
rownames(meta) <- cell_ids

# Short labels used downstream
meta$condition <- factor(c("normal" = "ND", "prediabetes syndrome" = "PD",
                           "type 2 diabetes mellitus" = "T2D")[as.character(meta$disease)],
                         levels = c("ND", "PD", "T2D"))
meta$donor     <- droplevels(meta$Islet)
meta$chemistry <- factor(sub("^10x 3' ", "", meta$assay))   # v2 / v3
meta$age       <- as.integer(sub("-year-old stage$", "", meta$development_stage))

## 4. Seurat object, log-normalised as in the GSE81608 pipeline
seu <- CreateSeuratObject(counts = counts, meta.data = meta,
                          project = "GSE221156", min.cells = 0, min.features = 0)
rm(counts)
seu <- NormalizeData(seu)

# Authors' UMAP from the beta-cell reintegration
umap <- t(h5read(h5ad_file, "obsm/X_umap"))
dimnames(umap) <- list(cell_ids, c("authorUMAP_1", "authorUMAP_2"))
seu[["umap_authors"]] <- CreateDimReducObject(umap, key = "authorUMAP_",
                                              assay = "RNA")
h5closeAll()

## 5. Summary
print(seu)
message("\nCells per condition:")
print(table(seu$condition))
message("\nDonors per condition:")
print(tapply(seu$donor, seu$condition, function(x) length(unique(x))))
message("\nBeta cells per donor:")
print(table(seu$donor, seu$condition)[, c("ND", "PD", "T2D")])
message("\nChemistry by condition:")
print(table(seu$chemistry, seu$condition))
message("\nAuthor beta subclusters by condition:")
print(table(seu$Clusters, seu$condition))

saveRDS(seu, rds_out)
message("\nSaved: ", rds_out)
