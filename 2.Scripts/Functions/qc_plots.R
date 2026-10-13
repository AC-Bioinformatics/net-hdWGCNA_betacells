### qc_plots(): per-group MAD filtering with kept/removed QC histograms
#
# Adapted from the APAP QC module (pipelines/qc_filtering/bin/qc_module.R):
#   1. For each group (e.g. donor / sample), work out two-sided MAD bounds
#      (median +/- n_mads*MAD) for every metric in `mad_filters`.
#   2. A cell is kept only if it passes every active filter.
#   3. Draw a histogram of every metric in `plot_metrics` on the *pre-filter*
#      data, colouring cells green (kept) / red (removed) by the FINAL
#      combined decision, with dashed vlines at the active bounds.
#   4. Write one PDF page per group, then a post-filter violin plot page.
#
# Usage:
#   source(file.path(Root, "2.Scripts", "Functions", "qc_plots.R"))
#   qc <- qc_plots(seu, group_by = "donor",
#                  mad_filters = list(nFeature_RNA = 5, percent.mt = 5),
#                  out_dir = Out_dir, log_fun = write_log)
#   seu      <- qc$seu        # all cells, with a qc_status column
#   seu_filt <- qc$seu_filt   # kept cells only
#   qc$bounds                 # data frame of bounds / pass counts per group

library(Seurat)
library(ggplot2)
library(patchwork)

`%||%` <- function(a, b) if (is.null(a)) b else a

# Two-sided MAD bounds: median +/- n_mads*MAD. Returns c(lower, upper).
# n_mads = NULL disables the bound (returns -Inf/Inf so it never filters anything).
mad_bounds <- function(x, n_mads, constant = 1.4826) {
  if (is.null(n_mads)) return(c(lower = -Inf, upper = Inf))
  med <- median(x, na.rm = TRUE)
  m   <- mad(x, na.rm = TRUE, constant = constant)
  c(lower = med - n_mads * m, upper = med + n_mads * m)
}

in_bounds <- function(x, bounds) {
  x > bounds["lower"] & x < bounds["upper"]
}

# Histogram of one QC variable, coloured by the FINAL keep/remove call,
# with optional dashed vlines for whichever bounds are active.
# vlines: named list, e.g. list(lower = list(x = 5, color = "steelblue"), ...)
qc_histogram <- function(df, var, vlines = list(), bins = 100,
                         title = var, xlab = var) {
  plot_df <- data.frame(value = df[[var]], status = df$qc_status)
  p <- ggplot(plot_df, aes(x = value, fill = status)) +
    geom_histogram(bins = bins, alpha = 0.85, color = NA, position = "stack") +
    scale_fill_manual(values = c(Kept = "#2ca02c", Removed = "#d62728"), name = "QC status") +
    theme_minimal() +
    labs(title = title, x = xlab, y = "Number of cells")

  line_labels <- c()
  for (nm in names(vlines)) {
    vl <- vlines[[nm]]
    if (!is.null(vl$x) && is.finite(vl$x)) {
      p <- p + geom_vline(xintercept = vl$x, color = vl$color,
                          linewidth = 0.9, linetype = "dashed")
      line_labels <- c(line_labels, paste0(nm, " = ", signif(vl$x, 4)))
    }
  }
  if (length(line_labels) > 0) {
    p <- p + labs(caption = paste(line_labels, collapse = "   |   "))
  }
  p
}

# seu          : Seurat object with the QC metrics already in meta.data
# group_by     : meta.data column to work out MAD bounds within (and one PDF page per value)
# mad_filters  : named list, metric = n_mads, e.g. list(nFeature_RNA = 5, percent.mt = 5)
#                Leave a metric out (or set it to NULL) to switch that filter off.
# plot_metrics : metrics to draw histograms (and violins) for; missing ones are skipped
# out_dir      : where the PDF goes
# pdf_name     : PDF file name (set to NULL to skip writing the PDF)
# log_fun      : function used for logging, e.g. write_log; defaults to message()
# bins, ncol, width, height : histogram bins, panels per row, PDF page size (inches)
# mad_constant : MAD-to-SD consistency constant (1.4826 = R's mad() default)
qc_plots <- function(seu,
                     group_by = "donor",
                     mad_filters = list(nFeature_RNA = 5, percent.mt = 5),
                     plot_metrics = c("nFeature_RNA", "nCount_RNA",
                                      "percent.mt", "percent.ribo", "percent.hb"),
                     out_dir = ".",
                     pdf_name = "QC_MAD_histograms.pdf",
                     log_fun = NULL,
                     bins = 100,
                     ncol = 3,
                     width = 17,
                     height = 11,
                     mad_constant = 1.4826) {
  log_fun <- log_fun %||% function(...) message(paste(...))
  mad_filters <- Filter(Negate(is.null), mad_filters)

  md <- seu[[]]
  if (!group_by %in% colnames(md)) {
    stop("group_by column '", group_by, "' not found in meta.data")
  }
  missing_filt <- setdiff(names(mad_filters), colnames(md))
  if (length(missing_filt) > 0) {
    stop("mad_filters metric(s) not in meta.data: ", paste(missing_filt, collapse = ", "))
  }
  missing_plot <- setdiff(plot_metrics, colnames(md))
  if (length(missing_plot) > 0) {
    warning("Skipping plot_metrics not in meta.data: ", paste(missing_plot, collapse = ", "))
    plot_metrics <- intersect(plot_metrics, colnames(md))
  }

  grp <- as.character(md[[group_by]])
  if (anyNA(grp)) {
    warning(sum(is.na(grp)), " cells have NA in '", group_by, "' and will be removed")
  }
  groups <- if (is.factor(md[[group_by]])) {
    intersect(levels(md[[group_by]]), grp)
  } else {
    sort(unique(na.omit(grp)))
  }

  log_fun("\n=====", format(Sys.time()), "- MAD QC filtering =====")
  log_fun("Grouped by:", group_by, "(", length(groups), "groups )")
  if (length(mad_filters) == 0) {
    log_fun("No MAD filters active - all cells kept")
  }
  for (metric in names(mad_filters)) {
    log_fun(" ", metric, "nmads:", mad_filters[[metric]])
  }

  keep_all <- rep(FALSE, nrow(md))
  all_hist_plots <- list()
  bounds_rows <- list()

  for (g in groups) {
    idx <- !is.na(grp) & grp == g
    md_g <- md[idx, , drop = FALSE]

    # Calculate bounds PER GROUP and combine all filters
    keep <- rep(TRUE, nrow(md_g))
    bounds <- list()
    log_fun("\n", group_by, ":", g, "| cells:", nrow(md_g))
    for (metric in names(mad_filters)) {
      b <- mad_bounds(md_g[[metric]], mad_filters[[metric]], constant = mad_constant)
      pass <- in_bounds(md_g[[metric]], b)
      pass[is.na(pass)] <- FALSE
      if (isTRUE(b["lower"] == b["upper"])) {
        warning("MAD of ", metric, " is 0 in ", group_by, " '", g,
                "' - every cell in this group fails that filter")
      }
      keep <- keep & pass
      bounds[[metric]] <- b
      log_fun("  ", metric, "bounds:", round(b["lower"], 2), "-", round(b["upper"], 2),
              "| pass:", sum(pass))
      bounds_rows[[length(bounds_rows) + 1]] <- data.frame(
        group = g, metric = metric, nmads = mad_filters[[metric]],
        lower = unname(b["lower"]), upper = unname(b["upper"]),
        pass = sum(pass), n_cells = nrow(md_g)
      )
    }
    md_g$qc_status <- ifelse(keep, "Kept", "Removed")
    keep_all[idx] <- keep
    log_fun("  Combined (final) pass:", sum(keep), "/", nrow(md_g),
            "(", round(100 * sum(keep) / nrow(md_g), 1), "%)")

    # Histograms for this group (pre-filter data, coloured by final call)
    p <- list()
    for (metric in plot_metrics) {
      b <- bounds[[metric]]
      vlines <- if (is.null(b)) list() else
        list(lower = list(x = b["lower"], color = "steelblue"),
             upper = list(x = b["upper"], color = "darkorange"))
      p[[metric]] <- qc_histogram(md_g, metric, vlines = vlines, bins = bins,
                                  title = paste(g, "-", metric), xlab = metric)
    }
    all_hist_plots[[g]] <- p
  }

  # Add final QC status to the object and filter
  seu$qc_status <- ifelse(keep_all, "Kept", "Removed")
  seu_filt <- subset(seu, cells = colnames(seu)[keep_all])

  log_fun("\nTotal cells before filtering:", ncol(seu))
  log_fun("Total cells after filtering :", ncol(seu_filt),
          "(", round(100 * ncol(seu_filt) / ncol(seu), 1), "%)")

  # Post-filter violin plot, same layout as the pre-filter one
  vln_title <- if (length(mad_filters) > 0) "Post-filter QC metrics by" else "QC metrics by"
  vln_post <- VlnPlot(seu_filt, features = plot_metrics, group.by = group_by,
                      pt.size = 0.1, ncol = ncol) +
    plot_annotation(title = paste(vln_title, group_by))

  # Histograms - one page per group, then the post-filter violin
  if (!is.null(pdf_name)) {
    pdf_path <- file.path(out_dir, pdf_name)
    pdf(pdf_path, width = width, height = height)
    for (g in groups) {
      print(wrap_plots(all_hist_plots[[g]], ncol = ncol) +
              plot_annotation(title = paste("QC Histograms -", g),
                              theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))))
    }
    print(vln_post)
    dev.off()
    log_fun("QC plots written to:", pdf_path)
  }

  invisible(list(
    seu = seu,
    seu_filt = seu_filt,
    bounds = do.call(rbind, bounds_rows),
    hist_plots = all_hist_plots,
    vln_post = vln_post
  ))
}


### ORA helpers for 6.module_pathways.R: plots and text summary for a clusterProfiler
### enricher() result. Requires clusterProfiler and enrichplot to be loaded.
#
# Usage:
#   source(file.path(Root, "2.Scripts", "Functions", "qc_plots.R"))
#   save_plot(p, file.path(plots_dir, "name_dotplot"))           # -> .png + .pdf
#   edo <- ora_plots(ora, analysis_name, plots_dir, fold_change = kme_vector)
#   write_ora_summary(txt_path, ora.res, params = list(...), ...)

# Write one plot as .png and .pdf (file_stem = path without extension). A plot that errors
# while being built or drawn (e.g. emapplot / treeplot with too few terms) is skipped: its
# device is closed and any partial files are removed.
save_plot <- function(expr, file_stem,
                      png_width = 3200, png_height = 2400, png_res = 300,
                      pdf_width = 12,   pdf_height = 9) {
  files <- paste0(file_stem, c(".png", ".pdf"))
  n_dev <- length(dev.list())
  ok <- tryCatch({
    p <- expr

    png(files[1], width = png_width, height = png_height, res = png_res, type = "cairo")
    print(p)
    dev.off()

    pdf(files[2], width = pdf_width, height = pdf_height)
    print(p)
    dev.off()
    TRUE
  }, error = function(e) {
    while (length(dev.list()) > n_dev) dev.off()
    unlink(files)
    message("Skipped: ", basename(file_stem), " (", conditionMessage(e), ")")
    FALSE
  })

  if (ok) message("Saved: ", basename(file_stem), "  (.png + .pdf)")
  invisible(ok)
}

# Standard ORA plots for an enrichResult: barplot, dotplot, emapplot, treeplot and cnetplot
# (genes coloured by fold_change, e.g. kME; NULL = uncoloured). subtitle is added to the
# barplot, dotplot and cnetplot. Returns the pairwise_termsim() object (NULL if it couldn't
# be built).
ora_plots <- function(ora, analysis_name, plots_dir, show_category = 20, fold_change = NULL,
                      subtitle = NULL) {
  stem <- function(x) file.path(plots_dir, paste0(analysis_name, "_", x))
  n_terms <- nrow(as.data.frame(ora))
  edo <- NULL

  # Suppress known harmless deprecation warnings from enrichplot
  withCallingHandlers(
    {
      save_plot(
        barplot(ora, showCategory = show_category) +
          labs(title = paste0(analysis_name, " - barplot"), subtitle = subtitle) +
          theme(axis.text.y = element_text(size = 6)),
        stem("barplot")
      )

      save_plot(
        dotplot(ora, showCategory = show_category) +
          labs(title = paste0(analysis_name, " - dotplot"), subtitle = subtitle) +
          theme(axis.text.y = element_text(size = 6)),
        stem("dotplot")
      )

      # Term-similarity plots need at least two terms
      if (n_terms >= 2) {
        edo <- tryCatch(pairwise_termsim(ora), error = function(e) NULL)
      }
      if (!is.null(edo)) {
        save_plot(emapplot(edo), stem("emapplot"))
        save_plot(treeplot(edo, cluster.params = list(n = min(5, n_terms))), stem("treeplot"),
                  png_width = 4200, pdf_width = 16)
      }

      save_plot(
        cnetplot(ora, showCategory = 5, foldChange = fold_change) +
          labs(title = paste0(analysis_name, " - gene-concept network",
                              if (!is.null(fold_change)) " (genes coloured by kME)"),
               subtitle = subtitle),
        stem("cnetplot")
      )
    },
    warning = function(w) {
      if (grepl("size|linewidth|by = .Count", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
  invisible(edo)
}

# Text summary of an ORA run.
# ora.res      : as.data.frame() of the significant enrichResult
# nominal.res  : optional as.data.frame() of the nominal-p terms, listed after the significant ones
# params       : named list printed in the PARAMETERS block, e.g. list(Module = "purple")
# n_background : universe size;  input_genes: genes tested (first top_genes listed)
# n_tested     : number of gene sets tested (all terms, before the cut-off)
# files        : named list of output paths printed at the end
write_ora_summary <- function(txt_path, ora.res, params, n_background, input_genes,
                              n_tested = NA, files = list(), top_genes = 10,
                              nominal.res = NULL, nominal_label = "nominal p < 0.05") {
  rule <- function(ch) cat(strrep(ch, 80), "\n", sep = "")
  print_terms <- function(df) {
    for (i in seq_len(nrow(df))) {
      cat(sprintf("[%2d] %s\n", i, df$ID[i]))
      cat(sprintf("     GeneRatio : %s\n",   df$GeneRatio[i]))
      cat(sprintf("     BgRatio   : %s\n",   df$BgRatio[i]))
      cat(sprintf("     p-value   : %.4e\n", df$pvalue[i]))
      cat(sprintf("     p.adjust  : %.4e\n", df$p.adjust[i]))
      cat(sprintf("     q.value   : %.4e\n", df$qvalue[i]))
      cat(sprintf("     Count     : %d\n",   df$Count[i]))
      cat(sprintf("     Genes     : %s\n\n", df$geneID[i]))
    }
  }

  sink(txt_path)
  on.exit(sink(), add = TRUE)
  rule("=")
  cat("ORA ANALYSIS SUMMARY\n")
  rule("=")
  cat(sprintf("%-16s: %s\n", "Date", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  rule("-")
  cat("PARAMETERS\n")
  rule("-")
  for (nm in names(params)) cat(sprintf("%-16s: %s\n", nm, paste(params[[nm]], collapse = ", ")))
  rule("-")
  cat("GENE LISTS\n")
  rule("-")
  cat(sprintf("%-16s: %d\n", "Background genes", n_background))
  cat(sprintf("%-16s: %d\n", "Input genes", length(input_genes)))
  cat(sprintf("%-16s: %s\n", paste("Top", top_genes), paste(head(input_genes, top_genes), collapse = ", ")))
  rule("-")
  cat("RESULTS OVERVIEW\n")
  rule("-")
  cat("Gene sets tested :", n_tested, "\n")
  cat("Enriched pathways:", nrow(ora.res), "\n")
  if (!is.null(nominal.res)) cat(sprintf("Terms with %s: %d\n", nominal_label, nrow(nominal.res)))

  if (nrow(ora.res) > 0) {
    cat("\n--- TOP 10 ENRICHED PATHWAYS (by p.adjust) ---\n\n")
    print_terms(head(ora.res[order(ora.res$p.adjust), ], 10))

    cat("--- FULL RESULTS TABLE ---\n\n")
    print(ora.res[, c("ID", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "Count")],
          row.names = FALSE)
  }

  if (!is.null(nominal.res) && nrow(nominal.res) > 0) {
    cat(sprintf("\n--- TOP 20 TERMS WITH %s (by p-value) ---\n\n", toupper(nominal_label)))
    print_terms(head(nominal.res[order(nominal.res$pvalue), ], 20))
  }

  cat("\n")
  rule("=")
  cat("Output files\n")
  rule("=")
  for (nm in names(files)) cat(sprintf("%-12s: %s\n", nm, files[[nm]]))
  invisible(txt_path)
}


### write_network_summary(): text summary of one module network (7.network_visualisation.R)
#
# Usage:
#   write_network_summary(txt_path, hub_tab, params = list(...), n_edges, n_pairs, files = list(...))
#
# hub_tab : one row per module gene with kME, kWithin, kWithin_scaled, kOut, kDiff, degree,
#           betweenness, hub and hallmark columns, ordered by kME
# params  : named list printed in the PARAMETERS block
# n_edges : edges kept in the thresholded network;  n_pairs: all within-module gene pairs
# files   : named list of output paths printed at the end
write_network_summary <- function(txt_path, hub_tab, params, n_edges, n_pairs,
                                  files = list(), top_n = 20) {
  rule <- function(ch) cat(strrep(ch, 80), "\n", sep = "")

  sink(txt_path)
  on.exit(sink(), add = TRUE)
  rule("=")
  cat("MODULE NETWORK SUMMARY\n")
  rule("=")
  cat(sprintf("%-20s: %s\n", "Date", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  rule("-")
  cat("PARAMETERS\n")
  rule("-")
  for (nm in names(params)) cat(sprintf("%-20s: %s\n", nm, paste(params[[nm]], collapse = ", ")))
  rule("-")
  cat("NETWORK\n")
  rule("-")
  cat(sprintf("%-20s: %d\n", "Genes", nrow(hub_tab)))
  cat(sprintf("%-20s: %d of %d (density %.3f)\n", "Edges kept", n_edges, n_pairs, n_edges / n_pairs))
  cat(sprintf("%-20s: %.1f (median %d)\n", "Mean degree", mean(hub_tab$degree), as.integer(median(hub_tab$degree))))
  cat(sprintf("%-20s: %d\n", "Genes without edges", sum(hub_tab$degree == 0)))
  cat(sprintf("%-20s: %.2f\n", "cor(kME, kWithin)", cor(hub_tab$kME, hub_tab$kWithin)))
  cat(sprintf("%-20s: %.2f (Spearman)\n", "cor(kME, degree)", cor(hub_tab$kME, hub_tab$degree, method = "spearman")))
  rule("-")
  cat(sprintf("TOP %d GENES BY kME\n", top_n))
  rule("-")
  cols <- c("gene_name", "kME", "kWithin_scaled", "kOut", "kDiff", "degree", "betweenness", "hub", "hallmark")
  top <- head(hub_tab[, cols], top_n)
  top[, c("kME", "kWithin_scaled", "kOut", "kDiff", "betweenness")] <-
    lapply(top[, c("kME", "kWithin_scaled", "kOut", "kDiff", "betweenness")], signif, 3)
  print(top, row.names = FALSE)
  rule("-")
  cat(sprintf("TOP %d GENES BY kWithin\n", 10))
  rule("-")
  print(head(hub_tab[order(-hub_tab$kWithin), c("gene_name", "kME", "kWithin_scaled", "degree", "hallmark")], 10),
        row.names = FALSE)
  rule("-")
  cat("GENES PER HALLMARK LABEL\n")
  rule("-")
  print(table(hub_tab$hallmark))
  cat("\n")
  rule("=")
  cat("Output files\n")
  rule("=")
  for (nm in names(files)) cat(sprintf("%-16s: %s\n", nm, files[[nm]]))
  invisible(txt_path)
}


### write_community_summary(): text summary of Leiden communities within a module
### (test.module_communities.R)
#
# res_stats : one row per resolution (n_communities_*, modularity_*, stability_ARI, sizes)
# nodes     : one row per gene with community, kME, kWithin_community, rank_in_community,
#             hub and hallmark columns
# params    : named list printed in the PARAMETERS block
# null_q    : modularity of the shuffled-weight null runs
# files     : named list of output paths printed at the end
write_community_summary <- function(txt_path, res_stats, nodes, params, null_q,
                                    files = list(), top_n = 10) {
  rule <- function(ch) cat(strrep(ch, 80), "\n", sep = "")

  sink(txt_path)
  on.exit(sink(), add = TRUE)
  rule("=")
  cat("MODULE COMMUNITY SUMMARY (LEIDEN)\n")
  rule("=")
  cat(sprintf("%-22s: %s\n", "Date", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  rule("-")
  cat("PARAMETERS\n")
  rule("-")
  for (nm in names(params)) cat(sprintf("%-22s: %s\n", nm, paste(params[[nm]], collapse = ", ")))
  rule("-")
  cat("RESOLUTIONS\n")
  rule("-")
  out <- as.data.frame(res_stats)
  num <- sapply(out, is.numeric)
  out[num] <- lapply(out[num], signif, 3)
  print(out, row.names = FALSE)
  cat(sprintf("\nNull modularity (shuffled weights): mean %.4f, max %.4f\n", mean(null_q), max(null_q)))
  rule("-")
  cat("COMMUNITIES AT THE CHOSEN RESOLUTION\n")
  rule("-")
  print(table(nodes$community))
  for (cm in sort(unique(nodes$community))) {
    d <- nodes[nodes$community == cm, ]
    cat(sprintf("\n--- %s: %d genes, mean kME %.2f, hubs: %s ---\n", cm, nrow(d), mean(d$kME),
                paste(d$gene_name[d$hub %in% "hub"], collapse = ", ")))
    cat("Top genes by within-community connectivity:",
        paste(head(d$gene_name[order(d$rank_in_community)], top_n), collapse = ", "), "\n")
    cat("Top genes by kME:", paste(head(d$gene_name[order(-d$kME)], top_n), collapse = ", "), "\n")
    if (!all(is.na(d$hallmark))) {
      tab <- sort(table(d$hallmark[d$hallmark != "none"]), decreasing = TRUE)
      cat("Hallmark labels:", if (length(tab)) paste(names(tab), tab, sep = " = ", collapse = "; ") else "none", "\n")
    }
  }
  cat("\n")
  rule("=")
  cat("Output files\n")
  rule("=")
  for (nm in names(files)) cat(sprintf("%-12s: %s\n", nm, files[[nm]]))
  invisible(txt_path)
}
