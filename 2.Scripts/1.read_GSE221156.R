### Script to read the GSE221156 beta-cell h5ad (CELLxGENE, Bandesh et al. 2026)
### into a Seurat object, keeping everything in the file as it is:
###   raw.X (raw UMI counts)            -> "counts" layer
###   X (authors' normalised values)    -> "data" layer (so no NormalizeData())
###   obs (all cell metadata)           -> meta.data, every column unchanged
###   var (all gene metadata)           -> RNA assay feature metadata, plus Ensembl IDs
###   obsm/X_umap (authors' UMAP)       -> "umap_authors" reduction
###   var$vst.variable                  -> VariableFeatures()
###   uns (title, citation, schema...)  -> Misc(seu, "h5ad_uns")
### The only changes are the ones Seurat forces on gene names (no underscores,
### unique names); the original symbol and Ensembl ID are kept for every gene.
### Short labels (condition, donor, chemistry, age) are added as extra columns.
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

## Read a matrix (raw/X or X) as genes x cells.
## Sparse CSR (cells x genes) is the same layout as CSC (genes x cells);
## dense arrays come back from rhdf5 already transposed to genes x cells.
read_matrix <- function(path) {
  if (paste0(path, "/data") %in% h5_names) {
    attrs <- h5readAttributes(h5ad_file, path)
    if (!identical(attrs[["encoding-type"]], "csr_matrix")) {
      stop(path, " is ", attrs[["encoding-type"]], "; expected csr_matrix")
    }
    # sparseMatrix() (not new("dgCMatrix")) because h5ad doesn't guarantee sorted indices
    sparseMatrix(i    = as.integer(h5read(h5ad_file, paste0(path, "/indices"))),
                 p    = as.integer(h5read(h5ad_file, paste0(path, "/indptr"))),
                 x    = as.numeric(h5read(h5ad_file, paste0(path, "/data"))),
                 dims = as.integer(rev(attrs$shape)),
                 index1 = FALSE)
  } else {
    as(Matrix(h5read(h5ad_file, path), sparse = TRUE), "CsparseMatrix")
  }
}

## 1. Raw counts
counts <- read_matrix("raw/X")
n_genes <- nrow(counts)
n_cells <- ncol(counts)
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

## 2b. Authors' normalised values (X), lined up to the raw genes
norm <- read_matrix("X")
stopifnot(ncol(norm) == n_cells)
gene_idx <- match(raw_ids, var_ids)
if (anyNA(gene_idx)) {
  stop(sum(is.na(gene_idx)), " genes in raw/X are not in X - can't use X as the data layer")
}
norm <- norm[gene_idx, ]
dimnames(norm) <- dimnames(counts)
if ("feature_is_filtered" %in% h5readAttributes(h5ad_file, "var")[["column-order"]]) {
  n_filt <- sum(as.logical(read_column("var", "feature_is_filtered"))[gene_idx])
  message("Genes the authors zeroed out in X (feature_is_filtered): ", n_filt)
}
# expm1(X) summing to ~10,000 per cell means log1p(counts per 10k), same as Seurat's LogNormalize
check_cells <- seq_len(min(100, n_cells))
message("Authors' X: median per-cell sum of expm1(X) = ",
        signif(median(Matrix::colSums(expm1(norm[, check_cells]))), 4),
        " (~10000 = log1p(CP10k))")

## 3. Cell metadata: every obs column, unchanged (incl. the authors' nCount_RNA /
##    nFeature_RNA / percent.mt - restored after CreateSeuratObject)
obs_cols <- h5readAttributes(h5ad_file, "obs")[["column-order"]]
meta <- as.data.frame(lapply(setNames(obs_cols, obs_cols),
                             function(col) read_column("obs", col)),
                      check.names = FALSE)
rownames(meta) <- cell_ids

# Short labels used downstream (added columns; the originals are kept as they are)
meta$condition <- factor(c("normal" = "ND", "prediabetes syndrome" = "PD",
                           "type 2 diabetes mellitus" = "T2D")[as.character(meta$disease)],
                         levels = c("ND", "PD", "T2D"))
meta$donor     <- droplevels(meta$Islet)
meta$chemistry <- factor(sub("^10x 3' ", "", meta$assay))   # v2 / v3
meta$age       <- as.integer(sub("-year-old stage$", "", meta$development_stage))

## 4. Seurat object: raw counts + the authors' normalised values as "data"
seu <- CreateSeuratObject(counts = counts, meta.data = meta,
                          project = "GSE221156", min.cells = 0, min.features = 0)
LayerData(seu, assay = "RNA", layer = "data") <- norm
# Setting a layer makes Seurat recompute nCount/nFeature (from that layer); put the
# authors' values back (computed on their full 36,601-gene data)
seu$nCount_RNA   <- meta$nCount_RNA
seu$nFeature_RNA <- meta$nFeature_RNA
rm(counts, norm)

# Gene metadata: every var column (lined up to the raw genes), plus the Ensembl ID
# and the original symbol, since Seurat's gene names may differ for a few genes
var_cols <- h5readAttributes(h5ad_file, "var")[["column-order"]]
var_meta <- as.data.frame(lapply(setNames(var_cols, var_cols),
                                 function(col) read_column("var", col)[gene_idx]),
                          check.names = FALSE)
var_meta <- cbind(ensembl_id = raw_ids, var_meta)
rownames(var_meta) <- gene_names
seu[["RNA"]] <- AddMetaData(seu[["RNA"]], metadata = var_meta)

# Authors' variable genes from their beta-cell reintegration (var$vst.variable)
VariableFeatures(seu) <- gene_names[as.logical(var_meta$vst.variable)]
message("Authors' variable features: ", length(VariableFeatures(seu)))

# File-level notes (title, citation, schema version, ...)
uns_names <- h5ls(h5ad_file)
uns_names <- uns_names$name[uns_names$group == "/uns"]
Misc(seu, "h5ad_uns") <- lapply(setNames(uns_names, uns_names),
                                function(u) as.vector(h5read(h5ad_file, paste0("uns/", u))))

# Authors' UMAP from the beta-cell reintegration
umap <- t(h5read(h5ad_file, "obsm/X_umap"))
dimnames(umap) <- list(cell_ids, c("authorUMAP_1", "authorUMAP_2"))
seu[["umap_authors"]] <- CreateDimReducObject(umap, key = "authorUMAP_",
                                              assay = "RNA")
h5closeAll()

## 5. Structure report: everything in the object, written to 3.Results
##    (overwritten each run, so it always describes the current object)
dir.create("3.Results", showWarnings = FALSE)
structure_file <- "3.Results/1.GSE221156_seurat_structure.txt"
cat("", file = structure_file)
out <- function(...) cat(paste0(...), "\n", file = structure_file, append = TRUE, sep = "")
out_obj <- function(x) cat(capture.output(x), file = structure_file, append = TRUE, sep = "\n")

# One line per column: type, then levels (with counts) or numeric range
describe <- function(x) {
  if (is.factor(x) || is.character(x) || is.logical(x)) {
    tb <- sort(table(x, useNA = "ifany"), decreasing = TRUE)
    lv <- paste0(names(tb), " (", tb, ")")
    sprintf("%s, %d values: %s%s", class(x)[1], length(tb),
            paste(head(lv, 4), collapse = "; "), if (length(lv) > 4) "; ..." else "")
  } else if (all(is.na(x))) {
    sprintf("%s, all NA", class(x)[1])
  } else {
    sprintf("%s: median %.4g, range %.4g - %.4g%s", class(x)[1], median(x, na.rm = TRUE),
            min(x, na.rm = TRUE), max(x, na.rm = TRUE),
            if (anyNA(x)) sprintf(", %d NA", sum(is.na(x))) else "")
  }
}

out("GSE221156 beta-cell Seurat object - structure")
out("Source: ", h5ad_file, "  |  created ", format(Sys.time()))
out("\n==================== OVERVIEW ====================")
out_obj(seu)

out("\n==================== LAYERS (RNA assay) ====================")
for (l in Layers(seu)) {
  m <- LayerData(seu, assay = "RNA", layer = l)
  out(sprintf("%-7s %s, %d genes x %d cells, %.1f%% non-zero, values %.4g - %.4g, all integers: %s",
              l, class(m)[1], nrow(m), ncol(m), 100 * length(m@x) / prod(dim(m)),
              min(m@x), max(m@x), all(m@x == round(m@x))))
}
out("counts = raw UMI counts (raw/X); data = authors' normalised values (X):")
out(sprintf("  median per-cell sum of expm1(data) = %.5g (~10000 = log1p(counts per 10k))",
            median(Matrix::colSums(expm1(LayerData(seu, layer = "data")[, check_cells])))))

out("\n==================== REDUCTIONS ====================")
for (r in Reductions(seu)) {
  out(sprintf("%-14s %d cells x %d dims (%s)", r, nrow(Embeddings(seu, r)), ncol(Embeddings(seu, r)),
              paste(colnames(Embeddings(seu, r)), collapse = ", ")))
}

out("\n==================== VARIABLE FEATURES ====================")
out(length(VariableFeatures(seu)), " genes (authors' vst.variable). First 30:")
out(paste(head(VariableFeatures(seu), 30), collapse = ", "))

out("\n==================== CELL METADATA (seu[[]]) ====================")
out(ncol(seu), " cells x ", ncol(seu[[]]), " columns")
for (col in colnames(seu[[]])) out(sprintf("  %-42s %s", col, describe(seu[[]][[col]])))

out("\n==================== GENE METADATA (seu[['RNA']][[]]) ====================")
out(nrow(seu), " genes x ", ncol(seu[["RNA"]][[]]), " columns")
for (col in colnames(seu[["RNA"]][[]])) out(sprintf("  %-28s %s", col, describe(seu[["RNA"]][[]][[col]])))
renamed <- which(as.character(var_meta$feature_name) != rownames(var_meta))
out(length(renamed), " genes have a Seurat name different from feature_name (underscores / duplicates)",
    if (length(renamed) > 0) paste0(", e.g.: ", paste(head(paste0(var_meta$feature_name[renamed], " -> ",
                                                                rownames(var_meta)[renamed]), 5), collapse = ", ")))

out("\n==================== FILE NOTES (Misc(seu, 'h5ad_uns')) ====================")
for (u in names(Misc(seu, "h5ad_uns"))) out(sprintf("  %-34s %s", u, paste(Misc(seu, "h5ad_uns")[[u]], collapse = "; ")))

out("\n==================== CELLS, DONORS AND CHEMISTRY ====================")
out("Cells per condition:")
out_obj(table(seu$condition))
out("\nDonors per condition:")
out_obj(tapply(seu$donor, seu$condition, function(x) length(unique(x))))
out("\nBeta cells per donor:")
out_obj(table(seu$donor, seu$condition)[, c("ND", "PD", "T2D")])
out("\nChemistry by condition:")
out_obj(table(seu$chemistry, seu$condition))
out("\nAuthor beta subclusters by condition:")
out_obj(table(seu$Clusters, seu$condition))
message("Structure report: ", structure_file)

saveRDS(seu, rds_out)
message("\nSaved: ", rds_out)
