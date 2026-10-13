# Pathway analysis (ORA, clusterProfiler) of every hdWGCNA module (grey excluded). Modules that
# track T2D in '5.hdWGCNA_module_correlations.R' (10.Module_condition_donor.csv /
# 10.Module_trait_cor_donor.csv) are run first:
#   purple  - lower in T2D (donor-level BH p = 0.002), HbA1c r = -0.68
#   magenta - higher in T2D (donor-level BH p = 0.014), HbA1c r = 0.52
# Module-cluster correlations are in 10.Module_cluster_cor.csv.
#
# Over-representation analysis: the genes assigned to a module are tested against every
# gene that went into the network (the universe, grey included).
# Gene sets come from msigdbr. "islet_t2d" is every C2 / C5 / Hallmark set whose name matches
# an islet, beta-cell, insulin or diabetes keyword.
# Tables hold the BH-significant terms (_results.csv) and the nominal p < 0.05 terms
# (_nominal.csv); plots show the nominal terms, coloured by BH-adjusted p.
# Plot and summary helpers are in Functions/qc_plots.R.

# Path configuration
Root <- "/home/alan.culligan/network"
Out_dir <- file.path(Root, "3.Results")
setwd(Root)

# Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(msigdbr)
})
source(file.path(Root, "2.Scripts", "Functions", "qc_plots.R"))

# Module assignment and kME of every network gene (GetModules(seu) after
# ModuleConnectivity, written by script 5)
kme <- read.csv(file.path(Out_dir, "10.Modules_kME.csv"))
universe <- kme$gene_name
message("Network genes (universe): ", length(universe))
print(table(kme$module))

run_module_ora <- function(
  module,
  pathway_db      = c("wikipathways", "reactome", "gobp", "hallmark", "islet_t2d"),
  analysis_name   = NULL,
  species         = "Homo sapiens",
  kme_min         = NULL,
  p_adjust_method = "BH",
  p_cutoff_ora    = 0.05,
  q_cutoff_ora    = 0.2,
  p_cutoff_plot   = 0.05,
  min_count_plot  = 2,
  min_gs_size     = 10,
  max_gs_size     = 500,
  show_category   = 20,
  output_dir      = file.path(Out_dir, "11.Module_pathways")
) {
  # ── Resolve pathway database ─────────────────────────────────────────────────
  pathway_db <- match.arg(pathway_db)

  # MSigDB category / subcategory for each database (msigdbr < 10 argument names);
  # pattern keeps only the gene sets whose name matches it
  db_params <- switch(pathway_db,
    wikipathways = list(
      category    = "C2",
      subcategory = "CP:WIKIPATHWAYS",
      label       = "WikiPathways"
    ),
    reactome = list(
      category    = "C2",
      subcategory = "CP:REACTOME",
      label       = "Reactome"
    ),
    gobp = list(
      category    = "C5",
      subcategory = "GO:BP",
      label       = "GO_BiologicalProcess"
    ),
    hallmark = list(
      category    = "H",
      subcategory = NULL,
      label       = "Hallmark"
    ),
    islet_t2d = list(
      category    = c("C2", "C5", "H"),
      subcategory = NULL,
      pattern     = paste0("DIABET|INSULIN|ISLET|BETA_CELL|PANCREATIC_BETA|MATURITY_ONSET|",
                           "RESPONSE_TO_GLUCOSE|GLUCOSE_STIMULATED|HNF1A|HNF4A|PDX1|NEUROG3|",
                           "UNFOLDED_PROTEIN|ENDOPLASMIC_RETICULUM_STRESS"),
      label       = "Islet_T2D"
    )
  )

  if (is.null(analysis_name)) {
    analysis_name <- paste0(module, "_ORA_", db_params$label)
  }

  # ── Output directories ───────────────────────────────────────────────────────
  plots_dir  <- file.path(output_dir, analysis_name, "plots")
  tables_dir <- file.path(output_dir, analysis_name, "tables")

  dir.create(plots_dir,  showWarnings = FALSE, recursive = TRUE)
  dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

  message("\n=== ORA | ", module, " | ", db_params$label, " | ", analysis_name, " ===")
  message("Output root : ", file.path(output_dir, analysis_name))

  # ── 1. Select module genes ───────────────────────────────────────────────────
  # Module genes = genes assigned to the module in GetModules(seu), with their kME to it.
  # kme_min keeps only the module core (e.g. 0.3); NULL keeps every assigned gene.
  kme_col <- paste0("kME_", module)
  if (!kme_col %in% colnames(kme)) stop("No module called '", module, "' in 10.Modules_kME.csv")

  module_genes <- kme %>%
    filter(module == !!module) %>%
    transmute(gene_name, kME = .data[[kme_col]]) %>%
    arrange(desc(kME))
  if (!is.null(kme_min)) module_genes <- filter(module_genes, kME >= kme_min)
  write.csv(module_genes, file.path(tables_dir, paste0(analysis_name, "_module_genes.csv")),
            row.names = FALSE)

  module.genes <- module_genes$gene_name
  bkgd.genes   <- universe
  message("Module genes: ", length(module.genes))
  message("Background  : ", length(bkgd.genes))
  if (length(module.genes) == 0) stop("No genes in module '", module, "' with the current kme_min.")

  # ── 2. Retrieve gene sets ────────────────────────────────────────────────────
  message("Fetching gene sets: ", db_params$label, "...")
  gs.sets <- bind_rows(lapply(db_params$category, function(cat) {
    msigdbr(species = species, category = cat, subcategory = db_params$subcategory)
  }))
  if (!is.null(db_params$pattern)) gs.sets <- filter(gs.sets, grepl(db_params$pattern, gs_name))
  gs.t2g <- gs.sets[, c("gs_name", "gene_symbol")]
  message("Gene sets retrieved: ", length(unique(gs.t2g$gs_name)))

  # ── 3. ORA ───────────────────────────────────────────────────────────────────
  # Every term is tested and saved (cut-offs = 1), then the p / q cut-offs are applied
  message("Running ORA...")
  ora_all <- clusterProfiler::enricher(
    gene          = module.genes,
    TERM2GENE     = gs.t2g,
    universe      = bkgd.genes,
    minGSSize     = min_gs_size,
    maxGSSize     = max_gs_size,
    pAdjustMethod = p_adjust_method,
    pvalueCutoff  = 1,
    qvalueCutoff  = 1
  )
  if (is.null(ora_all)) {
    warning(analysis_name, ": no gene set overlaps the module genes; skipping.")
    return(invisible(list(results = data.frame(), nominal = data.frame(), module_genes = module_genes)))
  }

  ora_all.res <- as.data.frame(ora_all)
  all_path <- file.path(tables_dir, paste0(analysis_name, "_all_terms.csv"))
  write.csv(ora_all.res, file = all_path, row.names = FALSE)

  ora <- dplyr::filter(ora_all, pvalue < p_cutoff_ora, p.adjust < p_cutoff_ora, qvalue < q_cutoff_ora)
  ora.res <- as.data.frame(ora)
  message("Significantly enriched pathways: ", nrow(ora.res), " of ", nrow(ora_all.res), " tested")

  csv_path <- file.path(tables_dir, paste0(analysis_name, "_results.csv"))
  write.csv(ora.res, file = csv_path, row.names = FALSE)
  message("Results saved to: ", csv_path)

  # Nominal terms (raw p-value, at least min_count_plot module genes) for the plots
  ora_nom <- dplyr::filter(ora_all, pvalue < p_cutoff_plot, Count >= min_count_plot)
  ora_nom.res <- as.data.frame(ora_nom)
  nom_path <- file.path(tables_dir, paste0(analysis_name, "_nominal.csv"))
  write.csv(ora_nom.res, file = nom_path, row.names = FALSE)
  message("Nominal terms (p < ", p_cutoff_plot, ", >= ", min_count_plot, " genes): ", nrow(ora_nom.res))

  # ── 4. Save text summary ─────────────────────────────────────────────────────
  txt_path <- file.path(tables_dir, paste0(analysis_name, "_summary.txt"))
  write_ora_summary(
    txt_path, ora.res,
    params = list(
      "Analysis name"    = analysis_name,
      "Module"           = module,
      "Pathway database" = db_params$label,
      "Species"          = species,
      "kME input"        = file.path(Out_dir, "10.Modules_kME.csv"),
      "kME minimum"      = if (is.null(kme_min)) "none (all assigned genes)" else kme_min,
      "Gene set size"    = paste(min_gs_size, "-", max_gs_size),
      "P-adjust method"  = p_adjust_method,
      "ORA p cutoff"     = p_cutoff_ora,
      "ORA q cutoff"     = q_cutoff_ora,
      "Plot p cutoff"    = paste(p_cutoff_plot, "(nominal), >=", min_count_plot, "genes")
    ),
    n_background = length(bkgd.genes),
    input_genes  = module.genes,
    n_tested     = nrow(ora_all.res),
    nominal.res  = ora_nom.res,
    files        = list("Tables dir" = tables_dir, "Plots dir" = plots_dir,
                        "CSV results" = csv_path, "Nominal" = nom_path,
                        "All terms" = all_path, "This file" = txt_path)
  )
  message("Saved: summary txt -> ", txt_path)

  # ── 5. Visualisations ────────────────────────────────────────────────────────
  edo <- NULL
  if (nrow(ora_nom.res) == 0) {
    warning(analysis_name, ": no nominal terms; skipping visualisations (see _all_terms.csv).")
  } else {
    message("Generating plots...")
    edo <- ora_plots(ora_nom, analysis_name, plots_dir, show_category = show_category,
                     fold_change = setNames(module_genes$kME, module_genes$gene_name),
                     subtitle = sprintf("Nominal p < %s, >= %d genes: %d terms (%d with BH p.adjust < %s); colour = p.adjust",
                                        p_cutoff_plot, min_count_plot, nrow(ora_nom.res),
                                        nrow(ora.res), p_cutoff_ora))
  }

  # ── 6. Return ────────────────────────────────────────────────────────────────
  invisible(list(
    ora          = ora,
    ora_all      = ora_all,
    ora_sim      = edo,
    results      = ora.res,
    nominal      = ora_nom.res,
    module_genes = module_genes,
    background   = bkgd.genes,
    pathway_db   = db_params$label,
    dirs         = list(plots = plots_dir, tables = tables_dir)
  ))
}

## Run: every module x each pathway database (purple and magenta first)
modules     <- union(c("purple", "magenta"), setdiff(names(sort(table(kme$module))), "grey"))
pathway_dbs <- c("wikipathways", "reactome", "gobp", "hallmark", "islet_t2d")

runs <- expand_grid(module = modules, pathway_db = pathway_dbs)
results <- pmap(runs, run_module_ora)
names(results) <- with(runs, paste(module, pathway_db, sep = "_"))

# Overview of nominal terms across all runs; significant = BH p.adjust < 0.05
overview <- bind_rows(lapply(seq_len(nrow(runs)), function(i) {
  df <- results[[i]]$nominal
  if (nrow(df) == 0) return(NULL)
  df %>%
    dplyr::select(ID, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count, geneID) %>%
    mutate(module = runs$module[i], pathway_db = runs$pathway_db[i],
           significant = p.adjust < 0.05, .before = 1)
}))
write.csv(overview, file.path(Out_dir, "11.Module_pathways_overview.csv"), row.names = FALSE)
message("\nSignificant (BH) / nominal terms per run:")
print(sapply(results, function(r) c(significant = nrow(r$results), nominal = nrow(r$nominal))))
message("Saved: ", file.path(Out_dir, "11.Module_pathways_overview.csv"))
