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
