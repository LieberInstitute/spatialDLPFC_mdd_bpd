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

JHPCE_HOST <- "jh"
JHPCE_REPO_ROOT <- "/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd"

DEG_FILE_SPECS <- list(
  smoothed_layer_adjusted = file.path(
    "processed-data", "07_dx_DE",
    "layer-adjusted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv"
  ),
  smoothed_layer_restricted = file.path(
    "processed-data", "07_dx_DE",
    "layer-restricted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv"
  ),
  seurat_layer_adjusted = file.path(
    "processed-data", "07_dx_DE",
    "layer-adjusted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv"
  ),
  seurat_layer_restricted = file.path(
    "processed-data", "07_dx_DE",
    "layer-restricted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv"
  )
)

AUTHOR_UNION_REL_PATH <- file.path(
  "raw-data", "SCENIC_aux", "tf_lists", "MBv_PRECAST-Seurat_F-test-adjp-05.txt"
)

SEURAT_CONTEXT_TO_DATASET_ID <- c(
  "Astro" = "astro",
  "Inhb" = "inhb",
  "L2.3" = "l2-3",
  "L4" = "l4",
  "L5" = "l5",
  "L6" = "l6",
  "Micro.Vasc" = "uvasc",
  "Oligo" = "oligo"
)

get_repo_root <- function(script_dir, args) {
  if (length(args) > 0) {
    normalizePath(args[[1]], mustWork = TRUE)
  } else {
    normalizePath(file.path(script_dir, "../.."), mustWork = TRUE)
  }
}

fetch_remote_file <- function(host, remote_path, local_path) {
  dir.create(dirname(local_path), recursive = TRUE, showWarnings = FALSE)

  test_status <- system2(
    "ssh",
    args = c("-o", "BatchMode=yes", host, "test", "-r", remote_path),
    stdout = FALSE,
    stderr = FALSE
  )
  if (!identical(test_status, 0L)) {
    stop("Remote DEG file is not readable via ", host, ": ", remote_path)
  }

  tmp_file <- tempfile(pattern = "deg-fetch-", tmpdir = dirname(local_path))
  err_file <- tempfile(pattern = "deg-fetch-err-", tmpdir = dirname(local_path))
  on.exit(unlink(c(tmp_file, err_file)), add = TRUE)

  status <- system2(
    "ssh",
    args = c("-o", "BatchMode=yes", host, "cat", remote_path),
    stdout = tmp_file,
    stderr = err_file
  )
  if (!identical(status, 0L)) {
    err <- character()
    if (file.exists(err_file)) {
      err <- readLines(err_file, warn = FALSE)
    }
    stop(
      "Failed to fetch DEG file from ", host, ": ", remote_path,
      if (length(err) > 0) paste0("\n", paste(err, collapse = "\n")) else ""
    )
  }
  if (!file.exists(tmp_file) || file.info(tmp_file)$size == 0) {
    stop("Fetched DEG file is empty: ", remote_path)
  }

  if (!file.rename(tmp_file, local_path)) {
    ok <- file.copy(tmp_file, local_path, overwrite = TRUE)
    if (!ok) stop("Failed to stage DEG file locally: ", local_path)
    unlink(tmp_file)
  }
}

ensure_local_file <- function(repo_root, rel_path, host = JHPCE_HOST, remote_root = JHPCE_REPO_ROOT) {
  local_path <- file.path(repo_root, rel_path)
  if (file.exists(local_path)) return(local_path)

  remote_path <- file.path(remote_root, rel_path)
  message("Staging missing file from ", host, ": ", rel_path)
  fetch_remote_file(host = host, remote_path = remote_path, local_path = local_path)
  local_path
}

load_deg_summary <- function(path, label) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  req <- c("gene_id", "gene_name", "adj.P.Val")
  missing_cols <- setdiff(req, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required DEG columns in ", label, " (", path, "): ",
      paste(missing_cols, collapse = ", ")
    )
  }

  df[df$adj.P.Val < 0.05, , drop = FALSE]
}

read_author_union_list <- function(path) {
  genes <- scan(path, what = "character", quiet = TRUE)
  sort(unique(genes[nzchar(genes)]))
}

collapse_gene_table <- function(df, label) {
  req <- c("gene_id", "gene_name")
  missing_cols <- setdiff(req, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing gene key columns in ", label, ": ", paste(missing_cols, collapse = ", "))
  }

  df2 <- df[, req, drop = FALSE]
  split_names <- split(df2$gene_name, df2$gene_id)
  bad_ids <- names(split_names)[vapply(split_names, function(x) length(unique(x)) != 1L, logical(1))]
  if (length(bad_ids) > 0) {
    stop(
      "Inconsistent gene_name values for gene_id(s) in ", label, ": ",
      paste(head(bad_ids, 10), collapse = ", ")
    )
  }

  df2[!duplicated(df2$gene_id), , drop = FALSE]
}

dataset_id_map_from_context_sets <- function(context_sets, context_to_dataset_id = SEURAT_CONTEXT_TO_DATASET_ID) {
  out <- list()
  for (ctx in names(context_to_dataset_id)) {
    ds_id <- unname(context_to_dataset_id[[ctx]])
    set_df <- context_sets[[ctx]]
    out[[ds_id]] <- set_df
    out[[paste0(ds_id, "_m")]] <- set_df
    out[[paste0(ds_id, "_f")]] <- set_df
  }
  out
}

validate_author_union <- function(deg_global, author_union_genes) {
  global_gene_names <- sort(unique(deg_global$gene_name))
  missing_from_author_union <- setdiff(global_gene_names, author_union_genes)
  extra_in_author_union <- setdiff(author_union_genes, global_gene_names)

  list(
    matches = length(missing_from_author_union) == 0L && length(extra_in_author_union) == 0L,
    gene_name_count = length(global_gene_names),
    author_gene_name_count = length(author_union_genes),
    missing_from_author_union = missing_from_author_union,
    extra_in_author_union = extra_in_author_union
  )
}

script_dir <- dirname(normalizePath(get_script_path(), mustWork = TRUE))
repo_root <- get_repo_root(script_dir, commandArgs(trailingOnly = TRUE))

deg_file_paths <- lapply(DEG_FILE_SPECS, function(rel_path) {
  ensure_local_file(repo_root = repo_root, rel_path = rel_path)
})
author_union_file <- ensure_local_file(repo_root = repo_root, rel_path = AUTHOR_UNION_REL_PATH)

deg_tables <- Map(
  f = function(path, label) load_deg_summary(path = path, label = label),
  path = deg_file_paths,
  label = names(deg_file_paths)
)
deg_tables <- deg_tables[names(DEG_FILE_SPECS)]

deg_global <- collapse_gene_table(
  do.call(
    rbind,
    lapply(unname(deg_tables), function(df) df[, c("gene_id", "gene_name"), drop = FALSE])
  ),
  label = "combined DEG summaries"
)

seurat_layer_adjusted_genes <- collapse_gene_table(
  deg_tables$seurat_layer_adjusted,
  label = "seurat_layer_adjusted"
)

deg_by_seurat_context <- lapply(names(SEURAT_CONTEXT_TO_DATASET_ID), function(context_name) {
  col_name <- paste0("n_ttest_sig_", context_name)
  lr_df <- deg_tables$seurat_layer_restricted
  if (!(col_name %in% names(lr_df))) {
    stop("Missing required context localization column in seurat_layer_restricted: ", col_name)
  }

  localized_lr <- lr_df[!is.na(lr_df[[col_name]]) & lr_df[[col_name]] > 0, , drop = FALSE]
  collapse_gene_table(
    rbind(seurat_layer_adjusted_genes, localized_lr[, c("gene_id", "gene_name"), drop = FALSE]),
    label = paste0("deg_by_seurat_context:", context_name)
  )
})
names(deg_by_seurat_context) <- names(SEURAT_CONTEXT_TO_DATASET_ID)

deg_by_dataset_id <- dataset_id_map_from_context_sets(deg_by_seurat_context)
author_union_genes <- read_author_union_list(author_union_file)
author_union_validation <- validate_author_union(deg_global, author_union_genes)

deg_validation <- list(
  source_counts = vapply(deg_tables, nrow, integer(1)),
  global_gene_id_count = nrow(deg_global),
  global_gene_name_count = length(unique(deg_global$gene_name)),
  seurat_lr_localized_gene_name_counts = vapply(
    names(SEURAT_CONTEXT_TO_DATASET_ID),
    function(context_name) {
      col_name <- paste0("n_ttest_sig_", context_name)
      sum(deg_tables$seurat_layer_restricted[[col_name]] > 0, na.rm = TRUE)
    },
    integer(1)
  ),
  seurat_context_gene_name_counts = vapply(
    deg_by_seurat_context,
    function(df) length(unique(df$gene_name)),
    integer(1)
  ),
  author_union = author_union_validation,
  files = c(
    deg_file_paths,
    author_union = author_union_file
  )
)

deg_ladj_file <- deg_file_paths$seurat_layer_adjusted
deg_lrestr_file <- deg_file_paths$seurat_layer_restricted
deg_ladj <- deg_tables$seurat_layer_adjusted
deg_lrestr <- deg_tables$seurat_layer_restricted

cat("Loaded author-recommended DEG summaries.\n")
cat("Repo root:", repo_root, "\n")
for (nm in names(deg_tables)) {
  cat(sprintf("  %s: %d rows\n", nm, nrow(deg_tables[[nm]])))
}
cat("Global union (gene_id):", deg_validation$global_gene_id_count, "\n")
cat("Global union (gene_name):", deg_validation$global_gene_name_count, "\n")
cat("Author union match:", deg_validation$author_union$matches, "\n")
cat("Per-context Seurat overlap set sizes (gene_name):\n")
for (ctx in names(deg_validation$seurat_context_gene_name_counts)) {
  ds_id <- unname(SEURAT_CONTEXT_TO_DATASET_ID[[ctx]])
  cat(sprintf("  %s (%s): %d\n", ctx, ds_id, deg_validation$seurat_context_gene_name_counts[[ctx]]))
}
