#!/usr/bin/env Rscript

get_script_path <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", ca, value = TRUE)
  if (length(m) > 0) return(sub("^--file=", "", m[[1]]))

  for (i in rev(seq_len(sys.nframe()))) {
    of <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(of)) return(of)
  }

  stop("Cannot determine script path")
}


## In the final DEG summary CSVs, adj.P.Val is the F-test BH-adjusted p-value.
## In the raw moderated t-test tables, adj.P.Val is also BH-adjusted, explicitly recomputed
##  globally in 03_plot-and-format_results.r and 03_plot-and-format_results.r

## So the DEG summary files already represent F-test-significant genes (omnibus test passed);
# additionally, we keep only genes with at least one significant moderated t-test contrast as well.
## i.e. genes that passed the omnibus F-test AND have at least one significant
## moderated t-test contrast.
## This yields 331 layer-adjusted rows and 203 layer-restricted rows in the current data.

## NOTE: the n_ttest_sig > 0 condition is added because 09_DEG_GRN/load_DEGs.r
## is a DEG loader that does the same filtering


load_deg_table <- function(path) {
  if (!file.exists(path)) stop("Missing DEG file: ", path)

  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  req <- c("adj.P.Val", "n_ttest_sig")
  missing_cols <- setdiff(req, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required DEG columns in ", path, ": ", paste(missing_cols, collapse = ", "))
  }

  df[df$adj.P.Val < 0.05 & df$n_ttest_sig > 0, , drop = FALSE]
  ## or could drop the n_ttest_sig filter
  # df[df$adj.P.Val < 0.05, , drop = FALSE]
}

script_dir <- dirname(normalizePath(get_script_path()))
repo_root <- if (length(commandArgs(trailingOnly = TRUE)) > 0) {
  normalizePath(commandArgs(trailingOnly = TRUE)[[1]])
} else {
  normalizePath(file.path(script_dir, "../.."))
}

deg_ladj_file <- file.path(
  repo_root,
  "processed-data",
  "07_dx_DE",
  "layer-adjusted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv"
)
deg_lrestr_file <- file.path(
  repo_root,
  "processed-data",
  "07_dx_DE",
  "layer-restricted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv"
)

deg_ladj <- load_deg_table(deg_ladj_file)
deg_lrestr <- load_deg_table(deg_lrestr_file)

cat("Loaded significant DEG summaries for seurat-pc30 strata.\n")
cat("Layer-adjusted:", nrow(deg_ladj), "rows\n")
cat("Layer-restricted:", nrow(deg_lrestr), "rows\n")
