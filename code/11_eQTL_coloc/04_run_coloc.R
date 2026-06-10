#!/usr/bin/env Rscript

## ---- setup ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(qs2)
  library(BiocParallel)
})

here::i_am(".git/HEAD")

repo_root <- here()
code_dir <- here("code", "11_eQTL_coloc")
source(file.path(code_dir, "utils.R"), chdir = FALSE)

## ---- config --------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

## set to "MDD" or "BPD" to restrict runs without CLI arguments.
RUN_DISORDER <- NULL

arg_value <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0L || hit[[1]] == length(args)) return(default)
  args[[hit[[1]] + 1L]]
}

arg_flag <- function(flag) flag %in% args

split_csv <- function(x, default) {
  if (is.null(x) || !nzchar(x)) return(default)
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

single_disorder_arg <- arg_value("--disorder")
if (!is.null(single_disorder_arg)) {
  single_disorder <- split_csv(single_disorder_arg, character())
  if (length(single_disorder) != 1L || !nzchar(single_disorder[[1]])) {
    stop("--disorder must specify exactly one disorder.")
  }
  if (!is.null(arg_value("--disorders"))) {
    message("--disorder provided; ignoring --disorders.")
  }
  disorders <- single_disorder
} else if (!is.null(RUN_DISORDER) && nzchar(RUN_DISORDER)) {
  disorders <- split_csv(RUN_DISORDER, character())
  if (length(disorders) != 1L || !nzchar(disorders[[1]])) {
    stop("RUN_DISORDER must be NULL or exactly one disorder.")
  }
} else {
  disorders <- split_csv(arg_value("--disorders"), c("MDD", "BPD"))
}
disorders <- vapply(disorders, gwas_check_disorder, character(1))
dataset_filter <- split_csv(arg_value("--datasets"), character())
chromosomes <- split_csv(arg_value("--chromosomes"), NULL)
dry_run <- arg_flag("--dry-run")
n_cores <- as.integer(arg_value("--n-cores", "4"))
if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L

tqtl_in_dir <- here("processed-data", "11_eQTL_coloc", "seurat", "tqtl_in")
tqtl_out_dir <- here("processed-data", "11_eQTL_coloc", "seurat", "tqtl_out")
coloc_dir <- here("processed-data", "11_eQTL_coloc", "seurat", "coloc")
plink2_prefix <- here("processed-data", "00_genotypes", "plink2", "merged_maf05")

## coloc GWAS slices are built on demand from the disorder-specific full BCF.
## each slice is dense over nominal eQTL variant positions and cached per
## disorder/domainCT. SI is the only GWAS row filter; no p-value filter is used.
cis_window <- 1000000L
si_min <- 0.8
min_snps <- 10L
min_abs_eqtl_z <- 2
sensitivity_rule <- "H4 > 0.8"
sensitivity_npoints <- 100L
priors <- list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)

required_packages <- c("arrow", "coloc", "BiocParallel", "qs2", "data.table")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "))
}

run_coloc_dir <- coloc_dir
if (!is.null(chromosomes)) {
  chr_tag <- gsub("[^A-Za-z0-9]+", "-", paste(chromosomes, collapse = "_"))
  run_coloc_dir <- file.path(coloc_dir, "validation", paste0("chromosomes_", chr_tag))
}
dir.create(run_coloc_dir, recursive = TRUE, showWarnings = FALSE)

## ---- manifest ------------------------------------------------------------
manifest <- load_eqtl_manifest(
  tqtl_in_dir = tqtl_in_dir,
  context_order = names(SEURAT_CONTEXT_TO_DATASET_ID),
  split_order = "all"
)

raw_manifest <- fread(file.path(tqtl_in_dir, "prep_manifest.csv"))
assert_cols(raw_manifest, c("dataset_id", "n_samples"), "prep_manifest.csv")
manifest <- merge(
  manifest,
  raw_manifest[, .(dataset_id, n_samples)],
  by = "dataset_id",
  all.x = TRUE,
  sort = FALSE
)
if (manifest[is.na(n_samples) | n_samples <= 0, .N] > 0L) {
  stop("Missing or invalid n_samples in coloc manifest.")
}

if (length(dataset_filter) > 0L) {
  manifest <- manifest[dataset_id %in% dataset_filter]
  if (nrow(manifest) == 0L) stop("No manifest rows matched --datasets.")
}
manifest[, context_order_rank := match(context, names(SEURAT_CONTEXT_TO_DATASET_ID))]
setorder(manifest, context_order_rank)
manifest[, context_order_rank := NULL]

message("Coloc datasets: ", paste(manifest$dataset_id, collapse = ", "))
message("Coloc disorders: ", paste(disorders, collapse = ", "))
message("Coloc output directory: ", run_coloc_dir)
if (!is.null(chromosomes)) message("Chromosome-restricted validation run: ", paste(chromosomes, collapse = ", "))
if (dry_run) {
  message("Dry run requested; no coloc outputs will be written.")
  quit(save = "no", status = 0)
}

## ---- run helpers ---------------------------------------------------------
result_path <- function(dis, dataset_id) {
  file.path(run_coloc_dir, dis, sprintf("coloc_%s.qs2", dataset_id))
}

meta_path <- function(dis, dataset_id) {
  file.path(run_coloc_dir, dis, sprintf("coloc_%s.runmeta.tsv.gz", dataset_id))
}

sensitivity_path <- function(dis, dataset_id) {
  file.path(run_coloc_dir, dis, sprintf("coloc_%s.sensitivity.tsv.gz", dataset_id))
}

complete_path <- function(dis, dataset_id) {
  file.path(run_coloc_dir, dis, sprintf("coloc_%s.complete", dataset_id))
}

validate_coloc_output_set <- function(out_path, runmeta_path, sens_path, dis, dataset_id) {
  dis_value <- dis
  dataset_id_value <- dataset_id
  needed <- c(result = out_path, metadata = runmeta_path, sensitivity = sens_path)
  missing <- needed[!file.exists(needed)]
  if (length(missing) > 0L) {
    stop("Missing coloc output file(s) for ", dis, "/", dataset_id, ": ", paste(names(missing), collapse = ", "))
  }
  bad_size <- needed[file.info(needed)$size <= 0]
  if (length(bad_size) > 0L) {
    stop("Empty coloc output file(s) for ", dis, "/", dataset_id, ": ", paste(names(bad_size), collapse = ", "))
  }

  meta <- fread(runmeta_path)
  assert_cols(meta, c("disorder", "dataset_id", "n_genes_saved"), basename(runmeta_path))
  if (nrow(meta) != 1L || !identical(as.character(meta$disorder[[1]]), dis) ||
      !identical(as.character(meta$dataset_id[[1]]), dataset_id)) {
    stop("Invalid coloc metadata for ", dis, "/", dataset_id)
  }

  res_list <- qs_read(out_path)
  if (!is.list(res_list) || length(res_list) != as.integer(meta$n_genes_saved[[1]])) {
    stop("Result object count does not match metadata for ", dis, "/", dataset_id)
  }

  sens <- fread(sens_path)
  assert_cols(sens, c("disorder", "dataset_id", "gene_id"), basename(sens_path))
  if (nrow(sens) > 0L && sens[disorder != dis_value | dataset_id != dataset_id_value, .N] > 0L) {
    stop("Sensitivity rows do not match metadata for ", dis, "/", dataset_id)
  }
  TRUE
}

write_complete_marker <- function(marker_path, out_path, runmeta_path, sens_path, dis, dataset_id) {
  validate_coloc_output_set(out_path, runmeta_path, sens_path, dis = dis, dataset_id = dataset_id)
  marker_tmp <- paste0(marker_path, ".tmp")
  writeLines(c(
    paste0("timestamp\t", as.character(Sys.time())),
    paste0("disorder\t", dis),
    paste0("dataset_id\t", dataset_id),
    paste0("result\t", out_path),
    paste0("metadata\t", runmeta_path),
    paste0("sensitivity\t", sens_path)
  ), marker_tmp)
  if (!file.rename(marker_tmp, marker_path)) {
    stop("Failed to write complete marker: ", marker_path)
  }
}

write_coloc_errors <- function(errors, err_path) {
  if (nrow(errors) > 0L) {
    fwrite(errors, err_path, sep = "\t", quote = FALSE, na = "NA")
    return(invisible(NULL))
  }
  con <- gzfile(err_path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines("disorder\tdataset_id\tgene_id\terror", con = con)
  invisible(NULL)
}

run_sensitivity_table <- function(result_list, dis, row) {
  dis_value <- dis
  dataset_id_value <- row$dataset_id
  context_value <- row$context
  split_value <- row$split
  if (length(result_list) == 0L) {
    return(data.table(
      disorder = character(),
      dataset_id = character(),
      context = character(),
      split = character(),
      gene_id = character(),
      cat = character(),
      PP3 = numeric(),
      PP4 = numeric(),
      PP34 = numeric(),
      sensitivity_error = character(),
      npoints = integer(),
      n_pass = integer(),
      frac_pass = numeric(),
      n_plausible = integer(),
      n_plausible_pass = integer(),
      frac_plausible_pass = numeric(),
      min_pass_p12 = numeric(),
      max_pass_p12 = numeric(),
      rule = character()
    ))
  }
  rows <- lapply(names(result_list), function(gene_id) {
    gene_id_value <- gene_id
    obj <- result_list[[gene_id]]
    sm <- as.list(obj$summary)
    pp3 <- as.numeric(sm[["PP.H3.abf"]])
    pp4 <- as.numeric(sm[["PP.H4.abf"]])
    sens <- summarize_coloc_sensitivity(
      obj,
      rule = sensitivity_rule,
      npoints = sensitivity_npoints
    )
    sens[, `:=`(
      disorder = dis_value,
      dataset_id = dataset_id_value,
      context = context_value,
      split = split_value,
      gene_id = gene_id_value,
      cat = coloc_assign_category_from_pp(pp3, pp4),
      PP3 = pp3,
      PP4 = pp4,
      PP34 = pp3 + pp4
    )]
    sens
  })
  out <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  front <- c("disorder", "dataset_id", "context", "split", "gene_id", "cat", "PP3", "PP4", "PP34")
  setcolorder(out, c(front, setdiff(names(out), front)))
  out[]
}

run_gene_workers <- function(genes, worker, dis, dataset_id) {
  named_genes <- setNames(genes, genes)
  if (n_cores <= 1L) return(lapply(named_genes, worker))

  ## coloc result objects are large; forked collection can fail if a worker dies.
  ## keep multicore optional and fall back to serial before writing any outputs.
  message(Sys.time(), " | ", dis, " | ", dataset_id, " | using ", n_cores, " coloc workers")
  bp <- BiocParallel::MulticoreParam(n_cores, stop.on.error = FALSE, progressbar = FALSE)
  res <- tryCatch(
    BiocParallel::bplapply(named_genes, worker, BPPARAM = bp),
    error = function(e) {
      message(
        Sys.time(), " | ", dis, " | ", dataset_id,
        " | multicore coloc failed; rerunning serial. Error: ", conditionMessage(e)
      )
      NULL
    }
  )
  if (is.null(res) || length(res) != length(named_genes)) return(lapply(named_genes, worker))
  res
}

run_coloc_dataset <- function(row, dis) {
  dataset_id <- row$dataset_id
  out_path <- result_path(dis, dataset_id)
  runmeta_path <- meta_path(dis, dataset_id)
  sens_path <- sensitivity_path(dis, dataset_id)
  done_path <- complete_path(dis, dataset_id)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

  if (file.exists(done_path) && file.exists(out_path) && file.exists(runmeta_path) && file.exists(sens_path)) {
    message(Sys.time(), " | ", dis, " | ", dataset_id, " | complete marker exists, skipping")
    return(invisible(NULL))
  }

  if (!file.exists(done_path) && file.exists(out_path) && file.exists(runmeta_path) && file.exists(sens_path)) {
    message(Sys.time(), " | ", dis, " | ", dataset_id, " | validating existing outputs before skip")
    write_complete_marker(done_path, out_path, runmeta_path, sens_path, dis = dis, dataset_id = dataset_id)
    return(invisible(NULL))
  }

  if (file.exists(out_path) && file.exists(runmeta_path) && !file.exists(sens_path)) {
    message(Sys.time(), " | ", dis, " | ", dataset_id, " | result exists; creating missing sensitivity table")
    res_list <- qs_read(out_path)
    sens <- run_sensitivity_table(res_list, dis = dis, row = row)
    fwrite(sens, sens_path, sep = "\t", quote = FALSE, na = "NA")
    write_complete_marker(done_path, out_path, runmeta_path, sens_path, dis = dis, dataset_id = dataset_id)
    return(invisible(NULL))
  }

  if (file.exists(out_path) || file.exists(runmeta_path)) {
    stop("Partial coloc output exists for ", dis, "/", dataset_id, ". Delete result and metadata before rerun.")
  }

  t0 <- Sys.time()
  message(Sys.time(), " | ", dis, " | ", dataset_id, " | reading nominal parquet")
  eqtl_dt <- read_coloc_nominal_dataset(
    dataset_id = dataset_id,
    tqtl_out_dir = tqtl_out_dir,
    cis_window = cis_window,
    chromosomes = chromosomes
  )

  ## prepare or reuse the disorder-specific GWAS cache for this domainCT.
  ## the queried regions are exactly the nominal eQTL variant positions that
  ## also exist in the PLINK2 genotype table, not pre-filtered GWAS hits.
  gwas_cache <- coloc_gwas_dataset_cache_file(
    dis = dis,
    dataset_id = dataset_id,
    coloc_dir = run_coloc_dir,
    si_min = si_min
  )
  message(Sys.time(), " | ", dis, " | ", dataset_id, " | extracting dense GWAS")
  gwas_dt <- extract_coloc_gwas_for_variants(
    dis = dis,
    variant_ids = unique(eqtl_dt$variant_id),
    plink2_prefix = plink2_prefix,
    out_file = gwas_cache,
    repo_root = repo_root,
    si_min = si_min
  )

  s_gwas <- coloc_case_fraction(gwas_dt)
  coloc_dt <- build_coloc_join_table(eqtl_dt, gwas_dt, plink2_prefix = plink2_prefix)
  genes <- sort(unique(as.character(coloc_dt$phenotype_id)))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  message(Sys.time(), " | ", dis, " | ", dataset_id, " | running coloc on ", length(genes), " gene loci")

  worker <- function(gene_id) {
    run_coloc_abf_one_gene_safe(
      gene_id = gene_id,
      coloc_dt = coloc_dt,
      n_eqtl = row$n_samples,
      s_gwas = s_gwas,
      p1 = priors$p1,
      p2 = priors$p2,
      p12 = priors$p12,
      min_snps = min_snps,
      min_abs_eqtl_z = min_abs_eqtl_z
    )
  }

  res_list <- run_gene_workers(genes = genes, worker = worker, dis = dis, dataset_id = dataset_id)

  err_idx <- which(vapply(res_list, inherits, logical(1), "coloc_err"))
  errors <- if (length(err_idx) > 0L) {
    data.table(
      disorder = dis,
      dataset_id = dataset_id,
      gene_id = names(res_list)[err_idx],
      error = vapply(res_list[err_idx], `[[`, character(1), ".error")
    )
  } else {
    data.table(disorder = character(), dataset_id = character(), gene_id = character(), error = character())
  }
  if (length(err_idx) > 0L) res_list <- res_list[-err_idx]

  null_idx <- vapply(res_list, is.null, logical(1))
  n_null <- sum(null_idx)
  res_list <- res_list[!null_idx]

  qs_save(res_list, out_path)
  sens <- run_sensitivity_table(res_list, dis = dis, row = row)
  fwrite(sens, sens_path, sep = "\t", quote = FALSE, na = "NA")

  err_path <- sub("\\.runmeta\\.tsv\\.gz$", ".errors.tsv.gz", runmeta_path)
  write_coloc_errors(errors, err_path)

  meta <- data.table(
    timestamp = as.character(Sys.time()),
    disorder = dis,
    dataset_id = dataset_id,
    context = row$context,
    split = row$split,
    cis_window = cis_window,
    chromosomes = if (is.null(chromosomes)) "" else paste(chromosomes, collapse = ","),
    gwas_si_min = si_min,
    gwas_cache = gwas_cache,
    out_path = out_path,
    sensitivity_path = sens_path,
    errors_path = err_path,
    n_eqtl_samples = row$n_samples,
    n_eqtl_rows = nrow(eqtl_dt),
    n_eqtl_genes = uniqueN(eqtl_dt$phenotype_id),
    n_eqtl_variants = uniqueN(eqtl_dt$variant_id),
    n_gwas_rows = nrow(gwas_dt),
    n_coloc_rows = nrow(coloc_dt),
    n_coloc_variants = uniqueN(coloc_dt$snp),
    n_genes_input = length(genes),
    n_genes_error = nrow(errors),
    n_genes_null = n_null,
    n_genes_saved = length(res_list),
    min_snps = min_snps,
    min_abs_eqtl_z = min_abs_eqtl_z,
    priors_p1 = priors$p1,
    priors_p2 = priors$p2,
    priors_p12 = priors$p12,
    case_fraction_s = s_gwas,
    elapsed_sec = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
  fwrite(meta, runmeta_path, sep = "\t", quote = FALSE, na = "NA")

  if (!file.exists(out_path) || !file.exists(runmeta_path) || !file.exists(sens_path)) {
    stop("Missing coloc output after run for ", dis, "/", dataset_id)
  }
  write_complete_marker(done_path, out_path, runmeta_path, sens_path, dis = dis, dataset_id = dataset_id)
  message(Sys.time(), " | ", dis, " | ", dataset_id, " | saved ", length(res_list), " coloc results")
  invisible(NULL)
}

## ---- execute -------------------------------------------------------------
for (dis in disorders) {
  for (i in seq_len(nrow(manifest))) {
    run_coloc_dataset(manifest[i], dis = dis)
  }
}

message("Coloc run complete.")
