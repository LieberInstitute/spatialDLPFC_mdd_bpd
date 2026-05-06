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

CUSTOM_CLUSTER_DEG_FILE_SPECS <- list(
  custom_layer_adjusted = file.path(
    "processed-data", "07_dx_DE",
    "layer-adjusted-pc3-age-nspots_custom-cluster_dx-sex_degs-F-test-t-test.csv"
  ),
  custom_layer_restricted = file.path(
    "processed-data", "07_dx_DE",
    "layer-restricted-pc3-age-nspots_custom-cluster_dx-sex_degs-F-test-t-test.csv"
  )
)

AUTHOR_UNION_REL_PATH <- file.path(
  "raw-data", "SCENIC_aux", "tf_lists", "MBv_PRECAST-Seurat_F-test-adjp-05.txt"
)

GENE_RANGES_REL_PATH <- file.path("processed-data", "ref", "granges.qs2")

DEFAULT_GENE_RANGES_SOURCE_RELS <- c(
  file.path(
    "processed-data", "06_pseudobulk", "Seurat",
    "spe_n119_pseudo-no-lowUMI_sample-seurat-pc30_norm-filt.Rdata"
  )
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

CUSTOM_CONTEXT_TO_DATASET_ID <- c(
  "Astro.L1" = "astro-l1",
  "Astro.Nrn" = "astro-nrn",
  "Inhb" = "inhb",
  "L2" = "l2",
  "L3" = "l3",
  "L4" = "l4",
  "L5" = "l5",
  "L6" = "l6",
  "Micro.Vasc" = "uvasc",
  "WM" = "wm"
)

DEG_SEX_PREFIX <- c(female = "F", male = "M")

DEG_TTEST_CONTRASTS <- c("NTC.MDD", "NTC.BPD", "MDD.BPD")

DEG_SPLITS <- c("all", "male", "female")

get_repo_root <- function(script_dir, args) {
  if (length(args) > 0) {
    normalizePath(args[[1]], mustWork = TRUE)
  } else {
    normalizePath(file.path(script_dir, "../.."), mustWork = TRUE)
  }
}

resolve_repo_root <- function(repo_root = NULL) {
  if (!is.null(repo_root)) return(normalizePath(repo_root, mustWork = TRUE))

  script_path <- tryCatch(get_script_path(), error = function(e) NULL)
  if (!is.null(script_path)) {
    script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
    return(get_repo_root(script_dir, commandArgs(trailingOnly = TRUE)))
  }

  if (file.exists(file.path(getwd(), "utils.R"))) {
    ## support source("utils.R") from this project directory.
    return(normalizePath(file.path(getwd(), "../.."), mustWork = TRUE))
  }

  stop("Cannot determine repo root; pass repo_root")
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

resolve_gene_range_source_file <- function(path, repo_root, host, remote_root) {
  if (grepl("^/", path)) {
    if (!file.exists(path)) stop("Missing gene range source file: ", path)
    return(normalizePath(path, mustWork = TRUE))
  }

  ensure_local_file(
    repo_root = repo_root,
    rel_path = path,
    host = host,
    remote_root = remote_root
  )
}

load_summarized_experiment_objects <- function(path) {
  obj_env <- new.env(parent = emptyenv())
  obj_names <- load(path, envir = obj_env)
  if (length(obj_names) == 0) stop("No objects found in ", path)

  is_se <- vapply(
    obj_names,
    function(x) inherits(get(x, envir = obj_env), "SummarizedExperiment"),
    logical(1)
  )
  if (!any(is_se)) {
    stop("No SummarizedExperiment object found in ", path)
  }

  lapply(obj_names[is_se], function(x) get(x, envir = obj_env))
}

validate_gene_ranges <- function(granges, label) {
  req <- c("gene_id", "gene_name")
  missing_cols <- setdiff(req, names(S4Vectors::mcols(granges)))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required gene range metadata in ", label, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  if (is.null(names(granges)) || anyNA(names(granges)) || any(!nzchar(names(granges)))) {
    stop("Missing GRanges names in ", label)
  }

  invisible(TRUE)
}

load_gene_ranges <- function(repo_root = NULL, source_files = NULL, granges_file = NULL,
                             verbose = TRUE, host = JHPCE_HOST,
                             remote_root = JHPCE_REPO_ROOT) {
  repo_root <- resolve_repo_root(repo_root)

  if (is.null(granges_file)) {
    granges_file <- file.path(repo_root, GENE_RANGES_REL_PATH)
  }
  granges_file <- normalizePath(granges_file, mustWork = FALSE)

  if (file.exists(granges_file)) {
    granges <- qs2::qs_read(granges_file)
    validate_gene_ranges(granges, granges_file)
    if (isTRUE(verbose)) cat("Loaded gene ranges:", granges_file, "\n")
    return(granges)
  }

  if (is.null(source_files)) source_files <- DEFAULT_GENE_RANGES_SOURCE_RELS
  source_paths <- vapply(
    source_files,
    resolve_gene_range_source_file,
    character(1),
    repo_root = repo_root,
    host = host,
    remote_root = remote_root
  )

  se_objects <- unlist(lapply(source_paths, load_summarized_experiment_objects), recursive = FALSE)
  all_granges <- lapply(se_objects, SummarizedExperiment::rowRanges)
  merged_granges <- do.call(c, unname(all_granges))

  ## preserve first occurrence when multiple source objects share gene IDs.
  merged_granges <- merged_granges[!duplicated(names(merged_granges))]
  validate_gene_ranges(merged_granges, paste(source_paths, collapse = ", "))

  dir.create(dirname(granges_file), recursive = TRUE, showWarnings = FALSE)
  qs2::qs_save(merged_granges, file = granges_file)
  if (isTRUE(verbose)) cat("Built gene ranges:", granges_file, "\n")

  merged_granges
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

ttest_sig_rows <- function(df, cols) {
  if (length(cols) == 0) return(rep(FALSE, nrow(df)))
  rowSums(df[, cols, drop = FALSE] == "padj<.05", na.rm = TRUE) > 0
}

validate_deg_view_columns <- function(seurat_lr_df,
                                      contexts = names(SEURAT_CONTEXT_TO_DATASET_ID),
                                      sex_prefix = DEG_SEX_PREFIX,
                                      contrasts = DEG_TTEST_CONTRASTS) {
  context_cols <- paste0("n_ttest_sig_", contexts)
  missing_context_cols <- setdiff(context_cols, names(seurat_lr_df))
  if (length(missing_context_cols) > 0) {
    stop("Missing Seurat context DEG columns: ", paste(missing_context_cols, collapse = ", "))
  }

  context_sex_cols <- unlist(lapply(contexts, function(context_name) {
    unlist(lapply(unname(sex_prefix), function(prefix) {
      paste0(context_name, "_", prefix, "_", contrasts, "_ttest")
    }), use.names = FALSE)
  }), use.names = FALSE)
  missing_context_sex_cols <- setdiff(context_sex_cols, names(seurat_lr_df))
  if (length(missing_context_sex_cols) > 0) {
    stop(
      "Missing Seurat context-by-sex DEG columns: ",
      paste(missing_context_sex_cols, collapse = ", ")
    )
  }

  invisible(TRUE)
}

make_deg_view_rows <- function(df, deg_view, context, split, deg_sex = "all") {
  df2 <- collapse_gene_table(df, label = paste(deg_view, context, split, deg_sex, sep = ":"))
  if (nrow(df2) == 0) {
    return(data.frame(
      deg_view = character(),
      context = character(),
      split = character(),
      deg_sex = character(),
      gene_id = character(),
      gene_name = character(),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    deg_view = deg_view,
    context = context,
    split = split,
    deg_sex = deg_sex,
    gene_id = df2$gene_id,
    gene_name = df2$gene_name,
    stringsAsFactors = FALSE
  )
}

sex_specific_deg_set <- function(deg_tables, sex, sex_prefix = DEG_SEX_PREFIX) {
  prefix <- unname(sex_prefix[[sex]])
  if (is.na(prefix)) stop("Unknown DEG sex: ", sex)

  out <- lapply(names(deg_tables), function(label) {
    df <- deg_tables[[label]]
    sex_cols <- grep(paste0("(^|_)", prefix, "_.*_ttest$"), names(df), value = TRUE)
    df[ttest_sig_rows(df, sex_cols), c("gene_id", "gene_name"), drop = FALSE]
  })

  collapse_gene_table(do.call(rbind, out), label = paste0("sex_specific:", sex))
}

build_deg_views <- function(deg_global, deg_tables,
                            contexts = names(SEURAT_CONTEXT_TO_DATASET_ID),
                            sex_prefix = DEG_SEX_PREFIX,
                            contrasts = DEG_TTEST_CONTRASTS,
                            splits = DEG_SPLITS) {
  seurat_lr_df <- deg_tables$seurat_layer_restricted
  validate_deg_view_columns(
    seurat_lr_df = seurat_lr_df,
    contexts = contexts,
    sex_prefix = sex_prefix,
    contrasts = contrasts
  )

  rows <- list()

  for (split_name in splits) {
    for (context_name in contexts) {
      rows[[length(rows) + 1L]] <- make_deg_view_rows(
        deg_global,
        deg_view = "broad_interaction",
        context = context_name,
        split = split_name
      )
    }
  }

  for (context_name in contexts) {
    col_name <- paste0("n_ttest_sig_", context_name)
    localized <- seurat_lr_df[
      !is.na(seurat_lr_df[[col_name]]) & seurat_lr_df[[col_name]] > 0,
      ,
      drop = FALSE
    ]
    for (split_name in splits) {
      rows[[length(rows) + 1L]] <- make_deg_view_rows(
        localized,
        deg_view = "context_localized",
        context = context_name,
        split = split_name
      )
    }
  }

  sex_sets <- lapply(names(sex_prefix), function(sex) {
    sex_specific_deg_set(deg_tables = deg_tables, sex = sex, sex_prefix = sex_prefix)
  })
  names(sex_sets) <- names(sex_prefix)

  for (context_name in contexts) {
    for (sex in names(sex_sets)) {
      for (split_name in c("all", sex)) {
        rows[[length(rows) + 1L]] <- make_deg_view_rows(
          sex_sets[[sex]],
          deg_view = "sex_specific",
          context = context_name,
          split = split_name,
          deg_sex = sex
        )
      }
    }
  }

  for (context_name in contexts) {
    for (sex in names(sex_prefix)) {
      prefix <- unname(sex_prefix[[sex]])
      sex_context_cols <- paste0(context_name, "_", prefix, "_", contrasts, "_ttest")
      sex_context <- seurat_lr_df[ttest_sig_rows(seurat_lr_df, sex_context_cols), , drop = FALSE]
      for (split_name in c("all", sex)) {
        rows[[length(rows) + 1L]] <- make_deg_view_rows(
          sex_context,
          deg_view = "context_and_sex_specific",
          context = context_name,
          split = split_name,
          deg_sex = sex
        )
      }
    }
  }

  long <- do.call(rbind, rows)
  rownames(long) <- NULL

  counts <- aggregate(
    gene_id ~ deg_view + context + split + deg_sex,
    data = long,
    FUN = function(x) length(unique(x)),
    na.action = NULL
  )
  names(counts)[names(counts) == "gene_id"] <- "n_genes"
  rownames(counts) <- NULL

  list(
    long = long,
    counts = counts,
    sex_sets = sex_sets
  )
}

build_custom_cluster_deg_views <- function(custom_by_context, broad_global, splits = c("all")) {
  rows <- list()
  for (split_name in splits) {
    for (context_name in names(custom_by_context)) {
      rows[[length(rows) + 1L]] <- make_deg_view_rows(
        custom_by_context[[context_name]],
        deg_view = "custom_primary",
        context = context_name,
        split = split_name
      )
      rows[[length(rows) + 1L]] <- make_deg_view_rows(
        broad_global,
        deg_view = "broad_PRECAST_Seurat_sensitivity",
        context = context_name,
        split = split_name
      )
    }
  }

  long <- do.call(rbind, rows)
  rownames(long) <- NULL

  counts <- aggregate(
    gene_id ~ deg_view + context + split + deg_sex,
    data = long,
    FUN = function(x) length(unique(x)),
    na.action = NULL
  )
  names(counts)[names(counts) == "gene_id"] <- "n_genes"
  rownames(counts) <- NULL

  list(
    long = long,
    counts = counts,
    sex_sets = list()
  )
}

load_standard_DEGs <- function(repo_root = NULL, verbose = TRUE,
                      host = JHPCE_HOST, remote_root = JHPCE_REPO_ROOT) {
  repo_root <- resolve_repo_root(repo_root)

  deg_file_paths <- lapply(DEG_FILE_SPECS, function(rel_path) {
    ensure_local_file(
      repo_root = repo_root,
      rel_path = rel_path,
      host = host,
      remote_root = remote_root
    )
  })
  author_union_file <- ensure_local_file(
    repo_root = repo_root,
    rel_path = AUTHOR_UNION_REL_PATH,
    host = host,
    remote_root = remote_root
  )

  deg_tables <- Map(
    f = function(path, label) load_deg_summary(path = path, label = label),
    path = deg_file_paths,
    label = names(deg_file_paths)
  )
  deg_tables <- deg_tables[names(DEG_FILE_SPECS)]

  ## global DEG support uses all four author-recommended F-test summaries.
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

  ## per-context eQTL overlap remains tied to current Seurat tensorQTL strata.
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
  deg_views <- build_deg_views(
    deg_global = deg_global,
    deg_tables = deg_tables,
    contexts = names(SEURAT_CONTEXT_TO_DATASET_ID)
  )
  author_union_genes <- read_author_union_list(author_union_file)

  ## author text list validates the reconstructed union; overlaps use tables above.
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
    deg_view_gene_counts = deg_views$counts,
    author_union = author_union_validation
  )

  files <- c(
    deg_file_paths,
    author_union = author_union_file
  )

  if (isTRUE(verbose)) {
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
  }

  list(
    tables = deg_tables,
    global = deg_global,
    by_seurat_context = deg_by_seurat_context,
    by_dataset_id = deg_by_dataset_id,
    views = deg_views,
    validation = deg_validation,
    files = files
  )
}

load_custom_cluster_DEGs <- function(repo_root = NULL, verbose = TRUE,
                                     host = JHPCE_HOST, remote_root = JHPCE_REPO_ROOT) {
  repo_root <- resolve_repo_root(repo_root)

  custom_file_paths <- lapply(CUSTOM_CLUSTER_DEG_FILE_SPECS, function(rel_path) {
    ensure_local_file(
      repo_root = repo_root,
      rel_path = rel_path,
      host = host,
      remote_root = remote_root
    )
  })
  custom_tables <- Map(
    f = function(path, label) load_deg_summary(path = path, label = label),
    path = custom_file_paths,
    label = names(custom_file_paths)
  )
  custom_tables <- custom_tables[names(CUSTOM_CLUSTER_DEG_FILE_SPECS)]

  broad_degs <- load_standard_DEGs(
    repo_root = repo_root,
    verbose = FALSE,
    host = host,
    remote_root = remote_root
  )

  custom_global <- collapse_gene_table(
    do.call(
      rbind,
      lapply(unname(custom_tables), function(df) df[, c("gene_id", "gene_name"), drop = FALSE])
    ),
    label = "custom-cluster DEG summaries"
  )

  custom_layer_adjusted_genes <- collapse_gene_table(
    custom_tables$custom_layer_adjusted,
    label = "custom_layer_adjusted"
  )

  custom_lr_df <- custom_tables$custom_layer_restricted
  custom_by_context <- lapply(names(CUSTOM_CONTEXT_TO_DATASET_ID), function(context_name) {
    col_name <- paste0("n_ttest_sig_", context_name)
    if (!(col_name %in% names(custom_lr_df))) {
      stop("Missing required custom-cluster DEG column: ", col_name)
    }

    localized_lr <- custom_lr_df[!is.na(custom_lr_df[[col_name]]) & custom_lr_df[[col_name]] > 0, , drop = FALSE]
    collapse_gene_table(
      rbind(custom_layer_adjusted_genes, localized_lr[, c("gene_id", "gene_name"), drop = FALSE]),
      label = paste0("deg_by_custom_context:", context_name)
    )
  })
  names(custom_by_context) <- names(CUSTOM_CONTEXT_TO_DATASET_ID)

  custom_by_dataset_id <- dataset_id_map_from_context_sets(
    custom_by_context,
    context_to_dataset_id = CUSTOM_CONTEXT_TO_DATASET_ID
  )
  custom_views <- build_custom_cluster_deg_views(
    custom_by_context = custom_by_context,
    broad_global = broad_degs$global,
    splits = c("all")
  )

  deg_validation <- list(
    source_counts = vapply(custom_tables, nrow, integer(1)),
    global_gene_id_count = nrow(custom_global),
    global_gene_name_count = length(unique(custom_global$gene_name)),
    broad_PRECAST_Seurat_gene_id_count = nrow(broad_degs$global),
    broad_PRECAST_Seurat_overlap_gene_id_count = length(intersect(custom_global$gene_id, broad_degs$global$gene_id)),
    custom_lr_localized_gene_name_counts = vapply(
      names(CUSTOM_CONTEXT_TO_DATASET_ID),
      function(context_name) {
        col_name <- paste0("n_ttest_sig_", context_name)
        sum(custom_lr_df[[col_name]] > 0, na.rm = TRUE)
      },
      integer(1)
    ),
    custom_context_gene_name_counts = vapply(
      custom_by_context,
      function(df) length(unique(df$gene_name)),
      integer(1)
    ),
    deg_view_gene_counts = custom_views$counts
  )

  files <- c(
    custom_file_paths,
    broad_sensitivity_files = broad_degs$files
  )

  if (isTRUE(verbose)) {
    cat("Loaded custom-cluster DEG summaries.\n")
    cat("Repo root:", repo_root, "\n")
    for (nm in names(custom_tables)) {
      cat(sprintf("  %s: %d rows\n", nm, nrow(custom_tables[[nm]])))
    }
    cat("Custom union (gene_id):", deg_validation$global_gene_id_count, "\n")
    cat("Broad PRECAST+Seurat union (gene_id):", deg_validation$broad_PRECAST_Seurat_gene_id_count, "\n")
    cat("Custom/broad overlap (gene_id):", deg_validation$broad_PRECAST_Seurat_overlap_gene_id_count, "\n")
    cat("Per-context custom overlap set sizes (gene_name):\n")
    for (ctx in names(deg_validation$custom_context_gene_name_counts)) {
      ds_id <- unname(CUSTOM_CONTEXT_TO_DATASET_ID[[ctx]])
      cat(sprintf("  %s (%s): %d\n", ctx, ds_id, deg_validation$custom_context_gene_name_counts[[ctx]]))
    }
  }

  list(
    tables = custom_tables,
    global = custom_global,
    by_custom_context = custom_by_context,
    by_dataset_id = custom_by_dataset_id,
    views = custom_views,
    validation = deg_validation,
    files = files,
    broad_sensitivity = broad_degs
  )
}

load_DEGs <- function(repo_root = NULL, mode = c("standard", "custom_cluster"), verbose = TRUE,
                      host = JHPCE_HOST, remote_root = JHPCE_REPO_ROOT) {
  mode <- match.arg(mode)
  if (identical(mode, "custom_cluster")) {
    return(load_custom_cluster_DEGs(
      repo_root = repo_root,
      verbose = verbose,
      host = host,
      remote_root = remote_root
    ))
  }

  load_standard_DEGs(
    repo_root = repo_root,
    verbose = verbose,
    host = host,
    remote_root = remote_root
  )
}

EQTL_REQUIRED_MAP_CIS_COLS <- c("phenotype_id", "variant_id", "qval")
EQTL_REQUIRED_INDEP_COLS <- c("phenotype_id", "variant_id", "rank", "pval_perm")

require_data_table <- function() {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("The data.table package is required")
  }
  invisible(TRUE)
}

path_relative_to <- function(path, root) {
  path_norm <- normalizePath(path, mustWork = FALSE)
  root_norm <- normalizePath(root, mustWork = TRUE)
  prefix <- paste0(root_norm, .Platform$file.sep)
  if (startsWith(path_norm, prefix)) return(sub(paste0("^", prefix), "", path_norm))
  path_norm
}

add_filename_suffix <- function(filename, suffix = "") {
  if (is.null(suffix) || !nzchar(suffix)) return(filename)
  sub("(\\.[^.]+)$", paste0(suffix, "\\1"), filename)
}

load_eqtl_manifest <- function(tqtl_in_dir, context_order = names(SEURAT_CONTEXT_TO_DATASET_ID),
                               split_order = DEG_SPLITS) {
  require_data_table()
  manifest_file <- file.path(tqtl_in_dir, "prep_manifest.csv")
  if (!file.exists(manifest_file)) stop("Missing tensorQTL manifest: ", manifest_file)

  manifest <- data.table::fread(manifest_file)
  req <- c("dataset_id", "seurat_label", "split", "status")
  missing_cols <- setdiff(req, names(manifest))
  if (length(missing_cols) > 0) {
    stop("Missing manifest columns in ", manifest_file, ": ", paste(missing_cols, collapse = ", "))
  }

  if (!"covariate_model" %in% names(manifest)) manifest[, covariate_model := "base"]
  manifest <- manifest[
    status == "prepared",
    .(dataset_id, context = seurat_label, split, covariate_model)
  ]
  manifest[, split := data.table::fifelse(
    split == "male",
    "male",
    data.table::fifelse(split == "female", "female", "all")
  )]

  if (anyDuplicated(manifest$dataset_id)) {
    stop("Duplicate dataset_id values in manifest: ", manifest_file)
  }
  unexpected_context <- setdiff(unique(manifest$context), context_order)
  if (length(unexpected_context) > 0) {
    stop("Unexpected manifest context(s): ", paste(unexpected_context, collapse = ", "))
  }
  unexpected_split <- setdiff(unique(manifest$split), split_order)
  if (length(unexpected_split) > 0) {
    stop("Unexpected manifest split(s): ", paste(unexpected_split, collapse = ", "))
  }

  manifest[order(match(split, split_order), match(context, context_order))]
}

read_eqtl_table <- function(dataset_id, tqtl_out_dir, suffix, required_cols,
                            missing_status, repo_root = NULL) {
  require_data_table()
  abs_path <- file.path(tqtl_out_dir, paste0(dataset_id, ".gene.", suffix))
  rel_path <- if (is.null(repo_root)) abs_path else path_relative_to(abs_path, repo_root)

  if (!file.exists(abs_path)) {
    return(list(status = missing_status, rel_path = NA_character_, dt = NULL))
  }
  if (file.info(abs_path)$size == 0) {
    stop("Empty tensorQTL result file: ", rel_path)
  }

  dt <- data.table::fread(abs_path)
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", rel_path, ": ", paste(missing_cols, collapse = ", "))
  }

  list(status = "ok", rel_path = rel_path, dt = dt)
}

read_map_cis <- function(dataset_id, tqtl_out_dir, repo_root = NULL,
                         required_cols = EQTL_REQUIRED_MAP_CIS_COLS) {
  read_eqtl_table(
    dataset_id = dataset_id,
    tqtl_out_dir = tqtl_out_dir,
    suffix = "map_cis.tab.gz",
    required_cols = required_cols,
    missing_status = "missing_map_cis",
    repo_root = repo_root
  )
}

read_map_independent <- function(dataset_id, tqtl_out_dir, repo_root = NULL,
                                 required_cols = EQTL_REQUIRED_INDEP_COLS) {
  read_eqtl_table(
    dataset_id = dataset_id,
    tqtl_out_dir = tqtl_out_dir,
    suffix = "map_independent.txt.gz",
    required_cols = required_cols,
    missing_status = "missing_map_independent",
    repo_root = repo_root
  )
}

collapse_egenes <- function(dt) {
  require_data_table()
  unique(dt[, .(gene_id = phenotype_id)])
}

collapse_gene_list <- function(x) {
  x <- sort(unique(x[!is.na(x) & nzchar(x)]))
  if (length(x) == 0) "" else paste(x, collapse = ", ")
}

tag_eqtl_results <- function(dt, dataset_id, context, split, deg_dt,
                             g2sym, gwas_set, gwas1e6_set) {
  require_data_table()
  out <- data.table::copy(dt)
  deg_gene_ids <- unique(data.table::as.data.table(deg_dt)$gene_id)

  out[, `:=`(
    dataset_id = dataset_id,
    context = context,
    split = split,
    gene_id = phenotype_id,
    gene_name = g2sym[phenotype_id],
    dge = as.integer(phenotype_id %in% deg_gene_ids),
    gwas = as.integer(variant_id %in% gwas_set),
    gwas_1e6 = as.integer(variant_id %in% gwas1e6_set)
  )]

  out[]
}

complete_summary_grid <- function(context_order = names(SEURAT_CONTEXT_TO_DATASET_ID),
                                  split_order = DEG_SPLITS) {
  require_data_table()
  data.table::CJ(split = split_order, context = context_order, unique = TRUE)
}

fill_summary_missing <- function(out, count_cols, text_cols) {
  for (col in count_cols) data.table::set(out, which(is.na(out[[col]])), col, 0L)
  for (col in text_cols) data.table::set(out, which(is.na(out[[col]])), col, "")
  out
}

summarize_tagged_eqtls <- function(dt, context_order = names(SEURAT_CONTEXT_TO_DATASET_ID),
                                   split_order = DEG_SPLITS) {
  require_data_table()
  if (nrow(dt) == 0) {
    out <- complete_summary_grid(context_order = context_order, split_order = split_order)
    out[, `:=`(
      n_eGenes = 0L,
      n_GWAS = 0L,
      n_GWAS_1e6 = 0L,
      n_DEG = 0L,
      n_DEG_GWAS = 0L,
      GWAS = "",
      DEGS = "",
      DEG_GWAS = ""
    )]
    return(out[])
  }

  obs <- dt[, .(
    n_eGenes = data.table::uniqueN(phenotype_id),
    n_GWAS = data.table::uniqueN(phenotype_id[gwas == 1]),
    n_GWAS_1e6 = data.table::uniqueN(phenotype_id[gwas_1e6 == 1]),
    n_DEG = data.table::uniqueN(gene_id[dge == 1]),
    n_DEG_GWAS = data.table::uniqueN(phenotype_id[dge == 1 & gwas == 1]),
    GWAS = collapse_gene_list(gene_name[gwas == 1]),
    DEGS = collapse_gene_list(gene_name[dge == 1]),
    DEG_GWAS = collapse_gene_list(gene_name[dge == 1 & gwas == 1])
  ), by = .(split, context)]

  out <- merge(
    complete_summary_grid(context_order = context_order, split_order = split_order),
    obs,
    by = c("split", "context"),
    all.x = TRUE
  )
  fill_summary_missing(
    out,
    count_cols = c("n_eGenes", "n_GWAS", "n_GWAS_1e6", "n_DEG", "n_DEG_GWAS"),
    text_cols = c("GWAS", "DEGS", "DEG_GWAS")
  )[order(match(split, split_order), match(context, context_order))]
}

summarize_independent_signals <- function(dt, context_order = names(SEURAT_CONTEXT_TO_DATASET_ID),
                                          split_order = DEG_SPLITS) {
  require_data_table()
  if (nrow(dt) == 0) {
    return(summarize_tagged_eqtls(dt, context_order = context_order, split_order = split_order))
  }

  obs <- dt[, .(
    n_eGenes = .N,
    n_GWAS = sum(gwas == 1),
    n_GWAS_1e6 = sum(gwas_1e6 == 1),
    n_DEG = data.table::uniqueN(gene_id[dge == 1]),
    n_DEG_GWAS = sum(dge == 1 & gwas == 1),
    GWAS = collapse_gene_list(gene_name[gwas == 1]),
    DEGS = collapse_gene_list(gene_name[dge == 1]),
    DEG_GWAS = collapse_gene_list(gene_name[dge == 1 & gwas == 1])
  ), by = .(split, context)]

  out <- merge(
    complete_summary_grid(context_order = context_order, split_order = split_order),
    obs,
    by = c("split", "context"),
    all.x = TRUE
  )
  fill_summary_missing(
    out,
    count_cols = c("n_eGenes", "n_GWAS", "n_GWAS_1e6", "n_DEG", "n_DEG_GWAS"),
    text_cols = c("GWAS", "DEGS", "DEG_GWAS")
  )[order(match(split, split_order), match(context, context_order))]
}

summarize_eqtl_deg_views <- function(eqtl_dt, manifest, degs,
                                     context_order = names(SEURAT_CONTEXT_TO_DATASET_ID),
                                     split_order = DEG_SPLITS,
                                     signal_count = FALSE) {
  require_data_table()
  view_dt <- data.table::as.data.table(degs$views$long)
  view_dt <- view_dt[, .(
    deg_view, context, split, deg_sex, gene_id,
    deg_gene_name = gene_name
  )]

  view_counts <- view_dt[, .(
    n_DEG_reference = data.table::uniqueN(gene_id)
  ), by = .(deg_view, context, split, deg_sex)]

  template <- merge(
    manifest[, .(dataset_id, context, split)],
    view_counts,
    by = c("context", "split"),
    all.x = FALSE,
    all.y = TRUE,
    allow.cartesian = TRUE
  )

  totals <- eqtl_dt[, .(
    n_eQTL_records = .N,
    n_eGenes = if (signal_count) .N else data.table::uniqueN(phenotype_id),
    n_GWAS = if (signal_count) sum(gwas == 1) else data.table::uniqueN(phenotype_id[gwas == 1]),
    n_GWAS_1e6 = if (signal_count) sum(gwas_1e6 == 1) else data.table::uniqueN(phenotype_id[gwas_1e6 == 1])
  ), by = .(dataset_id, context, split)]

  overlap_dt <- merge(
    eqtl_dt,
    view_dt,
    by = c("context", "split", "gene_id"),
    all = FALSE,
    allow.cartesian = TRUE
  )

  overlaps <- overlap_dt[, .(
    n_eGene_DEG_overlap = data.table::uniqueN(gene_id),
    n_DEG_GWAS = if (signal_count) sum(gwas == 1) else data.table::uniqueN(gene_id[gwas == 1]),
    eGene_DEG_overlap_gene_ids = collapse_gene_list(gene_id),
    eGene_DEG_overlap_gene_names = collapse_gene_list(deg_gene_name),
    DEG_GWAS = collapse_gene_list(deg_gene_name[gwas == 1])
  ), by = .(dataset_id, context, split, deg_view, deg_sex)]

  out <- merge(template, totals, by = c("dataset_id", "context", "split"), all.x = TRUE)
  out <- merge(
    out,
    overlaps,
    by = c("dataset_id", "context", "split", "deg_view", "deg_sex"),
    all.x = TRUE
  )

  fill_summary_missing(
    out,
    count_cols = c(
      "n_eQTL_records", "n_eGenes", "n_GWAS", "n_GWAS_1e6",
      "n_eGene_DEG_overlap", "n_DEG_GWAS"
    ),
    text_cols = c("eGene_DEG_overlap_gene_ids", "eGene_DEG_overlap_gene_names", "DEG_GWAS")
  )

  view_order <- c(
    "custom_primary",
    "broad_PRECAST_Seurat_sensitivity",
    "broad_interaction",
    "context_localized",
    "sex_specific",
    "context_and_sex_specific"
  )
  out[order(
    match(split, split_order),
    match(context, context_order),
    match(deg_view, view_order),
    deg_sex
  )]
}

summarize_dataset_counts <- function(dataset_id, context, split, file_info) {
  require_data_table()
  if (identical(file_info$status, "missing_map_independent")) {
    return(data.table::data.table(
      dataset_id = dataset_id,
      context = context,
      split = split,
      status = file_info$status,
      map_independent_file = NA_character_,
      n_independent_signals = NA_integer_,
      n_eGenes = NA_integer_,
      max_rank = NA_integer_,
      lead_variant_count = NA_integer_
    ))
  }

  dt <- file_info$dt
  data.table::data.table(
    dataset_id = dataset_id,
    context = context,
    split = split,
    status = file_info$status,
    map_independent_file = file_info$rel_path,
    n_independent_signals = nrow(dt),
    n_eGenes = data.table::uniqueN(dt$phenotype_id),
    max_rank = as.integer(max(dt$rank, na.rm = TRUE)),
    lead_variant_count = data.table::uniqueN(dt$variant_id)
  )
}

summarize_dataset_overlap <- function(dataset_id, context, split, file_info, deg_dt) {
  require_data_table()
  n_context_degs <- nrow(deg_dt)

  if (identical(file_info$status, "missing_map_independent")) {
    return(data.table::data.table(
      dataset_id = dataset_id,
      context = context,
      split = split,
      status = file_info$status,
      map_independent_file = NA_character_,
      n_independent_signals = NA_integer_,
      n_eGenes = NA_integer_,
      n_context_DEGs = n_context_degs,
      n_eGene_DEG_overlap = NA_integer_,
      eGene_DEG_overlap_gene_ids = NA_character_,
      eGene_DEG_overlap_gene_names = NA_character_
    ))
  }

  dt <- file_info$dt
  egenes <- collapse_egenes(dt)
  overlap_dt <- merge(
    egenes,
    data.table::as.data.table(deg_dt),
    by = "gene_id",
    all = FALSE,
    sort = FALSE
  )

  data.table::data.table(
    dataset_id = dataset_id,
    context = context,
    split = split,
    status = file_info$status,
    map_independent_file = file_info$rel_path,
    n_independent_signals = nrow(dt),
    n_eGenes = nrow(egenes),
    n_context_DEGs = n_context_degs,
    n_eGene_DEG_overlap = data.table::uniqueN(overlap_dt$gene_id),
    eGene_DEG_overlap_gene_ids = collapse_gene_list(overlap_dt$gene_id),
    eGene_DEG_overlap_gene_names = collapse_gene_list(overlap_dt$gene_name)
  )
}

summarize_eqtl_analysis <- function(config, degs, g2sym, gwas_set, gwas1e6_set,
                                    repo_root,
                                    context_order = names(SEURAT_CONTEXT_TO_DATASET_ID),
                                    split_order = DEG_SPLITS) {
  require_data_table()
  manifest <- load_eqtl_manifest(
    tqtl_in_dir = config$tqtl_in_dir,
    context_order = context_order,
    split_order = split_order
  )
  if (!all(manifest$dataset_id %in% names(degs$by_dataset_id))) {
    missing_ids <- setdiff(manifest$dataset_id, names(degs$by_dataset_id))
    stop("Missing DEG context sets for dataset IDs: ", paste(missing_ids, collapse = ", "))
  }

  summary_rows <- lapply(seq_len(nrow(manifest)), function(i) {
    row <- manifest[i]
    dataset_id <- row$dataset_id[[1]]
    context <- row$context[[1]]
    split <- row$split[[1]]

    map_cis_info <- read_map_cis(dataset_id, tqtl_out_dir = config$tqtl_out_dir, repo_root = repo_root)
    indep_info <- read_map_independent(dataset_id, tqtl_out_dir = config$tqtl_out_dir, repo_root = repo_root)
    deg_dt <- data.table::as.data.table(degs$by_dataset_id[[dataset_id]])

    list(
      counts = summarize_dataset_counts(dataset_id, context, split, indep_info),
      overlap = summarize_dataset_overlap(dataset_id, context, split, indep_info, deg_dt),
      map_cis = if (identical(map_cis_info$status, "ok")) {
        tag_eqtl_results(map_cis_info$dt, dataset_id, context, split, deg_dt, g2sym, gwas_set, gwas1e6_set)
      } else {
        NULL
      },
      independent = if (identical(indep_info$status, "ok")) {
        tag_eqtl_results(indep_info$dt, dataset_id, context, split, deg_dt, g2sym, gwas_set, gwas1e6_set)
      } else {
        NULL
      }
    )
  })

  counts_dt <- data.table::rbindlist(lapply(summary_rows, `[[`, "counts"), use.names = TRUE)
  overlap_dt <- data.table::rbindlist(lapply(summary_rows, `[[`, "overlap"), use.names = TRUE)
  map_cis_all <- data.table::rbindlist(lapply(summary_rows, `[[`, "map_cis"), use.names = TRUE, fill = TRUE)
  indep_all <- data.table::rbindlist(lapply(summary_rows, `[[`, "independent"), use.names = TRUE, fill = TRUE)

  counts_dt <- counts_dt[order(match(split, split_order), match(context, context_order))]
  overlap_dt <- overlap_dt[order(match(split, split_order), match(context, context_order))]

  map_cis_summary <- summarize_tagged_eqtls(
    map_cis_all[qval < 0.05],
    context_order = context_order,
    split_order = split_order
  )

  parent_q <- map_cis_all[, .(phenotype_id, dataset_id, qval_parent = qval)]
  indep_f <- merge(indep_all, parent_q, by = c("phenotype_id", "dataset_id"), all = FALSE)
  indep_f <- indep_f[qval_parent < 0.05 & pval_perm < 0.05]
  map_independent_summary <- summarize_independent_signals(
    indep_f,
    context_order = context_order,
    split_order = split_order
  )

  map_cis_deg_view_summary <- summarize_eqtl_deg_views(
    map_cis_all[qval < 0.05],
    manifest = manifest,
    degs = degs,
    context_order = context_order,
    split_order = split_order
  )
  map_independent_deg_view_summary <- summarize_eqtl_deg_views(
    indep_f,
    manifest = manifest,
    degs = degs,
    context_order = context_order,
    split_order = split_order,
    signal_count = TRUE
  )
  map_independent_deg_view_overlap <- summarize_eqtl_deg_views(
    indep_all,
    manifest = manifest,
    degs = degs,
    context_order = context_order,
    split_order = split_order
  )

  split_to_suffix <- c(all = "all", male = "male", female = "female")
  counts_tables <- lapply(names(split_to_suffix), function(split_name) counts_dt[split == split_name])
  names(counts_tables) <- names(split_to_suffix)
  overlap_tables <- lapply(names(split_to_suffix), function(split_name) overlap_dt[split == split_name])
  names(overlap_tables) <- names(split_to_suffix)

  list(
    analysis = config$analysis,
    manifest = manifest,
    counts_dt = counts_dt,
    overlap_dt = overlap_dt,
    map_cis_all = map_cis_all,
    indep_all = indep_all,
    indep_f = indep_f,
    map_cis_summary = map_cis_summary,
    map_independent_summary = map_independent_summary,
    map_cis_deg_view_summary = map_cis_deg_view_summary,
    map_independent_deg_view_summary = map_independent_deg_view_summary,
    map_independent_deg_view_overlap = map_independent_deg_view_overlap,
    counts_tables = counts_tables,
    overlap_tables = overlap_tables
  )
}
