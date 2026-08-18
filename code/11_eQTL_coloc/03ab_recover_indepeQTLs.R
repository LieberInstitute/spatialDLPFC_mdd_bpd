#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(dplyr)
  library(here)
  library(pgenlibr)
})

here::i_am(".git/HEAD")

repo_root <- here()
code_dir <- here("code", "11_eQTL_coloc")
source(file.path(code_dir, "utils.R"), chdir = FALSE)

## recover variants that are statistically indistinguishable from the single
## index variant emitted for each lead or conditionally independent eQTL signal.
## outputs use separate filenames unless --promote-current-pairs explicitly
## replaces the canonical pairs table after preserving its one-index input.

parse_cli <- function(args) {
  out <- list(
    datasets = NULL,
    chromosomes = NULL,
    output_dir = NULL,
    current_pairs = NULL,
    independent = NULL,
    coloc_gated = NULL,
    promote_current_pairs = FALSE,
    dosage_r_tolerance = 1e-12
  )
  for (arg in args) {
    if (startsWith(arg, "--datasets=")) {
      out$datasets <- strsplit(sub("^--datasets=", "", arg), ",", fixed = TRUE)[[1]]
    } else if (startsWith(arg, "--chromosomes=")) {
      out$chromosomes <- strsplit(sub("^--chromosomes=", "", arg), ",", fixed = TRUE)[[1]]
    } else if (startsWith(arg, "--output-dir=")) {
      out$output_dir <- sub("^--output-dir=", "", arg)
    } else if (startsWith(arg, "--current-pairs=")) {
      out$current_pairs <- sub("^--current-pairs=", "", arg)
    } else if (startsWith(arg, "--independent=")) {
      out$independent <- sub("^--independent=", "", arg)
    } else if (startsWith(arg, "--coloc-gated=")) {
      out$coloc_gated <- sub("^--coloc-gated=", "", arg)
    } else if (identical(arg, "--promote-current-pairs")) {
      out$promote_current_pairs <- TRUE
    } else if (startsWith(arg, "--dosage-r-tolerance=")) {
      out$dosage_r_tolerance <- as.numeric(sub("^--dosage-r-tolerance=", "", arg))
    } else {
      stop("Unknown argument: ", arg)
    }
  }
  if (!is.null(out$chromosomes)) {
    out$chromosomes <- unique(ifelse(
      startsWith(out$chromosomes, "chr"),
      out$chromosomes,
      paste0("chr", out$chromosomes)
    ))
  }
  if (!is.finite(out$dosage_r_tolerance) || out$dosage_r_tolerance <= 0) {
    stop("--dosage-r-tolerance must be a positive finite number")
  }
  out
}

opts <- parse_cli(commandArgs(trailingOnly = TRUE))

tqtl_in_dir <- here("processed-data", "11_eQTL_coloc", "seurat", "tqtl_in")
tqtl_out_dir <- here("processed-data", "11_eQTL_coloc", "seurat", "tqtl_out")
tables_dir <- here("processed-data", "11_eQTL_coloc", "seurat", "tables")
coloc_tables_dir <- here("processed-data", "11_eQTL_coloc", "seurat", "coloc", "tables")
plink2_prefix <- here("processed-data", "00_genotypes", "plink2", "merged_maf05")
output_dir <- if (is.null(opts$output_dir)) tables_dir else opts$output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_paths <- list(
  current_pairs = if (is.null(opts$current_pairs)) {
    file.path(tables_dir, "map_significant_pairs.csv.gz")
  } else {
    opts$current_pairs
  },
  independent = if (is.null(opts$independent)) {
    file.path(tables_dir, "map_independent_significant.csv.gz")
  } else {
    opts$independent
  },
  manifest = file.path(tqtl_in_dir, "prep_manifest.csv"),
  pgen = paste0(plink2_prefix, ".pgen"),
  pvar = paste0(plink2_prefix, ".pvar"),
  psam = paste0(plink2_prefix, ".psam"),
  coloc_gated = if (is.null(opts$coloc_gated)) {
    file.path(coloc_tables_dir, "coloc_abf_results_gated.tsv.gz")
  } else {
    opts$coloc_gated
  }
)

missing_inputs <- names(input_paths)[!file.exists(unlist(input_paths))]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing required input(s): ",
    paste(sprintf("%s=%s", missing_inputs, unlist(input_paths)[missing_inputs]), collapse = ", ")
  )
}

assert_required_cols <- function(dt, cols, label) {
  missing <- setdiff(cols, names(dt))
  if (length(missing) > 0L) {
    stop("Missing required columns in ", label, ": ", paste(missing, collapse = ", "))
  }
}

variant_chr <- function(x) {
  sub(":.*$", "", as.character(x))
}

nominal_file <- function(dataset_id, chromosome) {
  file.path(
    tqtl_out_dir,
    sprintf("%s.gene.cis_qtl_pairs.%s.parquet", dataset_id, chromosome)
  )
}

read_nominal_gene_subset <- function(path, gene_ids) {
  cols <- c(
    "phenotype_id", "variant_id", "start_distance", "af", "ma_samples",
    "ma_count", "pval_nominal", "slope", "slope_se"
  )
  ds <- arrow::open_dataset(path, format = "parquet")
  out <- ds |>
    dplyr::filter(phenotype_id %in% gene_ids) |>
    dplyr::select(dplyr::all_of(cols)) |>
    dplyr::collect()
  as.data.table(out)
}

mean_impute_and_standardize <- function(x) {
  x <- as.matrix(x)
  for (j in seq_len(ncol(x))) {
    missing <- is.na(x[, j])
    if (any(missing)) {
      observed <- x[!missing, j]
      if (length(observed) == 0L) stop("A candidate genotype vector is entirely missing")
      x[missing, j] <- mean(observed)
    }
  }
  x <- sweep(x, 2L, colMeans(x), "-")
  norms <- sqrt(colSums(x * x))
  if (any(!is.finite(norms) | norms <= 0)) stop("A candidate genotype vector has zero variance")
  sweep(x, 2L, norms, "/")
}

dataset_expression_samples <- function(dataset_id) {
  path <- file.path(tqtl_in_dir, paste0(dataset_id, ".gene.expr.bed.gz"))
  if (!file.exists(path)) stop("Missing tensorQTL expression BED: ", path)
  header <- names(fread(path, nrows = 0L))
  if (length(header) < 5L) stop("Expression BED has no sample columns: ", path)
  header[5:length(header)]
}

front_order <- function(dt, cols) {
  setcolorder(dt, c(intersect(cols, names(dt)), setdiff(names(dt), cols)))
  dt[]
}

current_pairs <- fread(input_paths$current_pairs)
independent <- fread(input_paths$independent)
manifest <- fread(input_paths$manifest)

recovery_metadata_cols <- c(
  "signal_id", "signal_kind", "signal_rank", "index_variant_id",
  "is_index_variant", "was_in_current_pairs", "recovery_status",
  "n_signal_members", "tie_orientation", "dosage_r", "dosage_r2",
  "signal_chr", "signal_rank_label", "n_samples", "slope_scale",
  "member_raw_start_distance", "member_raw_af", "member_raw_ma_samples",
  "member_raw_ma_count", "member_raw_pval_nominal", "member_raw_slope",
  "member_raw_slope_se", "index_raw_start_distance",
  "index_raw_pval_nominal", "index_raw_slope", "index_raw_slope_se"
)
if (all(c("signal_id", "is_index_variant", "recovery_status") %in% names(current_pairs))) {
  message("Current-pairs input is already tie-expanded; using its reported-index rows as the recovery baseline")
  current_pairs <- current_pairs[is_index_variant == TRUE]
  current_pairs[, (intersect(recovery_metadata_cols, names(current_pairs))) := NULL]
}

pair_keys <- c("dataset_id", "context", "split", "gene_id", "variant_id")
assert_required_cols(
  current_pairs,
  c(pair_keys, "gene_name", "phenotype_id", "cis_supported", "indep_supported"),
  basename(input_paths$current_pairs)
)
assert_required_cols(independent, c(pair_keys, "rank"), basename(input_paths$independent))
assert_required_cols(manifest, c("dataset_id", "seurat_label", "split", "n_samples", "status"), basename(input_paths$manifest))

if (anyDuplicated(current_pairs[, ..pair_keys])) {
  stop("Current significant-pair input is not unique by ", paste(pair_keys, collapse = ", "))
}

rank_lut <- unique(independent[, c(pair_keys, "rank"), with = FALSE])
setnames(rank_lut, "rank", "signal_rank")
signals <- merge(current_pairs, rank_lut, by = pair_keys, all.x = TRUE, sort = FALSE)

bad_independent_ranks <- signals[indep_supported == TRUE & is.na(signal_rank)]
if (nrow(bad_independent_ranks) > 0L) {
  stop("Failed to recover an independent-signal rank for ", nrow(bad_independent_ranks), " current pair(s)")
}

signals[, `:=`(
  signal_kind = fifelse(indep_supported == TRUE, "conditionally_independent", "lead_cis_only"),
  signal_chr = variant_chr(variant_id)
)]
signals[, signal_rank_label := fifelse(
  signal_kind == "conditionally_independent",
  paste0("rank", signal_rank),
  "lead"
)]
signals[, signal_id := paste(dataset_id, gene_id, signal_rank_label, sep = "|")]

if (anyDuplicated(signals$signal_id)) {
  dup <- unique(signals[duplicated(signal_id) | duplicated(signal_id, fromLast = TRUE), signal_id])
  stop("Signal IDs are not unique: ", paste(head(dup, 20L), collapse = ", "))
}

manifest_all <- manifest[status == "prepared" & split == "all"]
manifest_all[, context := seurat_label]
manifest_all <- manifest_all[, .(dataset_id, context, split, n_samples)]
signals <- merge(signals, manifest_all, by = c("dataset_id", "context", "split"), all.x = TRUE, sort = FALSE)
if (anyNA(signals$n_samples)) stop("Missing manifest sample counts for one or more signals")

psam <- fread(input_paths$psam)
psam_id_col <- if ("#IID" %in% names(psam)) "#IID" else if ("IID" %in% names(psam)) "IID" else NA_character_
if (is.na(psam_id_col)) stop("Cannot find IID column in ", input_paths$psam)
psam_samples <- as.character(psam[[psam_id_col]])

if (!is.null(opts$datasets)) {
  unknown <- setdiff(opts$datasets, unique(signals$dataset_id))
  if (length(unknown) > 0L) stop("Unknown selected dataset(s): ", paste(unknown, collapse = ", "))
  signals <- signals[dataset_id %in% opts$datasets]
}
if (!is.null(opts$chromosomes)) {
  signals <- signals[signal_chr %in% opts$chromosomes]
}
if (nrow(signals) == 0L) stop("No signals remain after dataset/chromosome filtering")

selected_datasets <- sort(unique(signals$dataset_id))
selected_chromosomes <- sort(unique(signals$signal_chr))
dataset_samples <- setNames(lapply(selected_datasets, dataset_expression_samples), selected_datasets)
for (dataset_id in selected_datasets) {
  ds_value <- dataset_id
  samples <- dataset_samples[[dataset_id]]
  missing <- setdiff(samples, psam_samples)
  if (length(missing) > 0L) {
    stop("Expression samples missing from PSAM for ", dataset_id, ": ", paste(head(missing, 20L), collapse = ", "))
  }
  expected_n <- unique(signals[dataset_id == ds_value, n_samples])
  if (length(expected_n) != 1L || length(samples) != expected_n) {
    stop(
      "Manifest/expression sample-count mismatch for ", dataset_id,
      ": manifest=", paste(expected_n, collapse = ","), ", expression BED=", length(samples)
    )
  }
}
message("Selected datasets: ", paste(selected_datasets, collapse = ", "))
message("Selected chromosomes: ", paste(selected_chromosomes, collapse = ", "))
message("Current index-variant rows/signals: ", nrow(signals))

required_nominal <- unique(signals[, .(
  path = nominal_file(dataset_id, signal_chr),
  dataset_id,
  chromosome = signal_chr
)])
missing_nominal <- required_nominal[!file.exists(path)]
if (nrow(missing_nominal) > 0L) {
  stop(
    "Missing ", nrow(missing_nominal), " required full nominal parquet file(s). First missing file(s): ",
    paste(head(missing_nominal$path, 10L), collapse = ", ")
  )
}

candidate_parts <- vector("list", nrow(required_nominal))
for (i in seq_len(nrow(required_nominal))) {
  spec <- required_nominal[i]
  signal_subset <- signals[dataset_id == spec$dataset_id & signal_chr == spec$chromosome]
  gene_ids <- unique(signal_subset$gene_id)
  message(
    sprintf(
      "[%d/%d] reading %s for %d signal(s) across %d gene(s)",
      i, nrow(required_nominal), basename(spec$path), nrow(signal_subset), length(gene_ids)
    )
  )
  nominal <- read_nominal_gene_subset(spec$path, gene_ids)
  if (nrow(nominal) == 0L) stop("No nominal rows recovered from ", spec$path)

  index_lut <- signal_subset[, .(
    signal_id,
    gene_id,
    index_variant_id = variant_id,
    signal_kind,
    signal_rank
  )]
  index_raw <- merge(
    index_lut,
    nominal,
    by.x = c("gene_id", "index_variant_id"),
    by.y = c("phenotype_id", "variant_id"),
    all.x = TRUE,
    sort = FALSE
  )
  if (anyNA(index_raw$pval_nominal)) {
    missing_index <- index_raw[is.na(pval_nominal), paste(gene_id, index_variant_id, sep = "/")]
    stop("Index variant missing from raw nominal data: ", paste(head(missing_index, 20L), collapse = ", "))
  }
  if (anyDuplicated(index_raw$signal_id)) stop("Duplicate raw nominal index rows found in ", spec$path)

  index_raw <- index_raw[, .(
    signal_id,
    gene_id,
    index_variant_id,
    signal_kind,
    signal_rank,
    index_raw_start_distance = start_distance,
    index_raw_pval_nominal = pval_nominal,
    index_raw_slope = slope,
    index_raw_slope_se = slope_se
  )]
  candidates <- merge(index_raw, nominal, by.x = "gene_id", by.y = "phenotype_id", allow.cartesian = TRUE, sort = FALSE)

  ## exact equality here reproduces tensorQTL's tied nominal test statistic.
  ## genotype equivalence is verified below before any row is recovered.
  candidates <- candidates[pval_nominal == index_raw_pval_nominal]
  candidates[, `:=`(
    dataset_id = spec$dataset_id,
    context = signal_subset$context[[1]],
    split = signal_subset$split[[1]],
    chromosome = spec$chromosome
  )]
  setnames(candidates, "variant_id", "member_variant_id")
  candidate_parts[[i]] <- candidates
  rm(nominal, candidates, index_raw)
  gc(verbose = FALSE)
}

candidates <- rbindlist(candidate_parts, use.names = TRUE, fill = TRUE)
if (nrow(candidates) == 0L) stop("No nominal tie candidates were found")
if (any(!signals$variant_id %in% candidates$member_variant_id)) {
  stop("At least one index variant is absent from the nominal tie candidates")
}
message("Nominal same-p candidate memberships: ", nrow(candidates))

candidate_variant_ids <- unique(c(candidates$index_variant_id, candidates$member_variant_id))
message("Loading genotype dosages for ", length(candidate_variant_ids), " unique candidate variant(s)")

pvar <- fread(input_paths$pvar, sep = "\t", skip = "#CHROM", select = 1:5)
setnames(pvar, 1:5, c("CHROM", "POS", "ID", "REF", "ALT"))
pvar[, var_idx := .I]
pvar_subset <- pvar[match(candidate_variant_ids, ID)]
if (anyNA(pvar_subset$var_idx)) {
  missing <- candidate_variant_ids[is.na(pvar_subset$var_idx)]
  stop("Candidate variant(s) missing from PVAR: ", paste(head(missing, 20L), collapse = ", "))
}

pgen <- NewPgen(input_paths$pgen)
dosage <- as.matrix(ReadList(pgen, pvar_subset$var_idx, meanimpute = FALSE))
rm(pgen)
colnames(dosage) <- pvar_subset$ID
rownames(dosage) <- psam_samples
dosage_col <- setNames(seq_len(ncol(dosage)), colnames(dosage))
candidates[, dosage_r := NA_real_]
for (dataset_id in selected_datasets) {
  ds_value <- dataset_id
  sample_rows <- match(dataset_samples[[dataset_id]], rownames(dosage))
  dosage_z <- mean_impute_and_standardize(dosage[sample_rows, , drop = FALSE])
  for (index_id in unique(candidates[dataset_id == ds_value, index_variant_id])) {
    rows <- which(candidates$dataset_id == ds_value & candidates$index_variant_id == index_id)
    member_cols <- unname(dosage_col[candidates$member_variant_id[rows]])
    index_col <- unname(dosage_col[[index_id]])
    candidates[rows, dosage_r := as.numeric(crossprod(
      dosage_z[, index_col],
      dosage_z[, member_cols, drop = FALSE]
    ))]
  }
  rm(dosage_z)
}
rm(dosage)
gc(verbose = FALSE)
candidates[dosage_r > 1 & dosage_r < 1 + opts$dosage_r_tolerance, dosage_r := 1]
candidates[dosage_r < -1 & dosage_r > -1 - opts$dosage_r_tolerance, dosage_r := -1]
candidates[, dosage_r2 := dosage_r * dosage_r]
candidates[, exact_dosage_tie := abs(abs(dosage_r) - 1) <= opts$dosage_r_tolerance]

rejected_same_p <- candidates[exact_dosage_tie == FALSE]
members <- candidates[exact_dosage_tie == TRUE]
if (nrow(members) == 0L) stop("No exact dosage-equivalent signal members were found")
if (any(!signals$signal_id %in% members$signal_id)) {
  stop("At least one signal lost its index variant during exact-dosage verification")
}

message("Verified exact-dosage signal memberships: ", nrow(members))
message("Same-p candidates rejected because |dosage r| < 1: ", nrow(rejected_same_p))

members[, `:=`(
  tie_orientation = fifelse(dosage_r >= 0, 1L, -1L),
  is_index_variant = member_variant_id == index_variant_id
)]
members[, slope_scale := slope / index_raw_slope]
members[!is.finite(slope_scale), slope_scale := as.numeric(tie_orientation)]
members[, n_signal_members := .N, by = signal_id]
members[, recovery_status := fifelse(is_index_variant, "reported_index", "recovered_exact_tie")]

current_pair_lut <- unique(signals[, .(dataset_id, context, split, gene_id, current_variant_id = variant_id)])
members[, was_in_current_pairs := FALSE]
members[current_pair_lut, on = .(
  dataset_id,
  context,
  split,
  gene_id,
  member_variant_id = current_variant_id
), was_in_current_pairs := TRUE]

member_fields <- c(
  "signal_id", "index_variant_id", "member_variant_id",
  "n_signal_members", "is_index_variant", "was_in_current_pairs", "recovery_status",
  "tie_orientation", "dosage_r", "dosage_r2", "slope_scale",
  "start_distance", "af", "ma_samples", "ma_count", "pval_nominal", "slope", "slope_se",
  "index_raw_start_distance", "index_raw_pval_nominal", "index_raw_slope", "index_raw_slope_se"
)
member_map <- members[, ..member_fields]
setnames(
  member_map,
  c("start_distance", "af", "ma_samples", "ma_count", "pval_nominal", "slope", "slope_se"),
  c(
    "member_raw_start_distance", "member_raw_af", "member_raw_ma_samples", "member_raw_ma_count",
    "member_raw_pval_nominal", "member_raw_slope", "member_raw_slope_se"
  )
)

expanded <- merge(signals, member_map, by = "signal_id", allow.cartesian = TRUE, sort = FALSE)
expanded[, index_end_distance := end_distance]
expanded[, variant_id := member_variant_id]
expanded[is_index_variant == FALSE, `:=`(
  start_distance = member_raw_start_distance,
  end_distance = index_end_distance + (member_raw_start_distance - index_raw_start_distance),
  af = member_raw_af,
  ma_samples = member_raw_ma_samples,
  ma_count = member_raw_ma_count,
  slope = slope * slope_scale,
  slope_se = slope_se * abs(slope_scale)
)]
expanded[, c("index_end_distance", "member_variant_id") := NULL]

## refresh exact-variant GWAS annotations for every recovered member while
## retaining the gene-level DEG and GWAS-gene annotations from its index row.
gwas_matched <- load_matched_gwas_by_disorder(
  DEFAULT_GWAS_OVERLAP_DISORDERS,
  plink2_prefix = plink2_prefix,
  repo_root = repo_root,
  si_min = GWAS_MATCH_SI_MIN
)
strict_gwas_sets <- strict_gwas_sets_from_matched(gwas_matched, GWAS_STRICT_P_THRESHOLD)
exploratory_gwas <- gwasx_matched_from_gwas(gwas_matched, GWAS_EXPLORATORY_P_THRESHOLDS)

for (dis in names(strict_gwas_sets)) {
  strict_variant_col <- paste0(dis, "_gwasVar_strict")
  gene_col <- paste0(dis, "_gwasGene")
  strict_combined_col <- paste0(dis, "_gwas_strict")
  expanded[is_index_variant == FALSE, (strict_variant_col) := as.integer(variant_id %in% strict_gwas_sets[[dis]])]
  if (gene_col %in% names(expanded)) {
    expanded[is_index_variant == FALSE, (strict_combined_col) := as.integer(
      get(strict_variant_col) == 1L | get(gene_col) == 1L
    )]
  }
}

for (dis in names(exploratory_gwas)) {
  exp_variant_col <- paste0(dis, "_gwasVar_exp")
  gene_col <- paste0(dis, "_gwasGene")
  exp_combined_col <- paste0(dis, "_gwas_exp")
  p_col <- paste0(dis, "_gwasP")
  beta_col <- paste0(dis, "_gwasBeta")
  beta_se_col <- paste0(dis, "_gwasBetaSE")
  stats <- exploratory_gwas[[dis]][, .(
    variant_id,
    gwas_p = p,
    gwas_beta = beta,
    gwas_beta_se = beta_se
  )]
  expanded[is_index_variant == FALSE, (exp_variant_col) := as.integer(variant_id %in% stats$variant_id)]
  expanded[is_index_variant == FALSE, c(p_col, beta_col, beta_se_col) := list(
    stats$gwas_p[match(variant_id, stats$variant_id)],
    stats$gwas_beta[match(variant_id, stats$variant_id)],
    stats$gwas_beta_se[match(variant_id, stats$variant_id)]
  )]
  if (gene_col %in% names(expanded)) {
    expanded[is_index_variant == FALSE, (exp_combined_col) := as.integer(
      get(exp_variant_col) == 1L | get(gene_col) == 1L
    )]
  }
}

## verify that refreshing annotations did not change any existing index flags.
flag_cols <- intersect(
  c(
    "MDD_gwasVar_strict", "MDD_gwasVar_exp", "BD_gwasVar_strict", "BD_gwasVar_exp",
    "SCZD_gwasVar_strict"
  ),
  names(signals)
)
index_check <- expanded[is_index_variant == TRUE, c(pair_keys, flag_cols), with = FALSE]
original_check <- signals[, c(pair_keys, flag_cols), with = FALSE]
check <- merge(index_check, original_check, by = pair_keys, suffixes = c("_new", "_old"), all = TRUE)
for (col in flag_cols) {
  if (any(check[[paste0(col, "_new")]] != check[[paste0(col, "_old")]], na.rm = TRUE)) {
    stop("Refreshed index-variant annotations disagree with the current table for ", col)
  }
}

metadata_front <- c(
  "signal_id", "signal_kind", "signal_rank", "index_variant_id", "variant_id",
  "is_index_variant", "was_in_current_pairs", "recovery_status", "n_signal_members",
  "tie_orientation", "dosage_r", "dosage_r2",
  "dataset_id", "context", "split", "gene_id", "gene_name", "phenotype_id"
)
expanded <- front_order(expanded, metadata_front)
setorder(expanded, dataset_id, gene_id, signal_rank, signal_id, -is_index_variant, variant_id)

if (anyDuplicated(expanded[, ..pair_keys])) {
  stop("Recovered alternative table is not unique by ", paste(pair_keys, collapse = ", "))
}

independent_expanded <- expanded[indep_supported == TRUE]

signal_member_counts <- expanded[, .(
  n_signal_members = .N,
  n_hidden_members = sum(was_in_current_pairs == FALSE),
  n_hidden_unique_variants = uniqueN(variant_id[was_in_current_pairs == FALSE])
), by = .(dataset_id, context, split, gene_id, gene_name, signal_id, signal_kind, signal_rank, index_variant_id)]

summary_by_context <- signal_member_counts[, .(
  n_current_signal_rows = .N,
  n_signals_with_exact_ties = sum(n_signal_members > 1L),
  n_genes_with_exact_ties = uniqueN(gene_id[n_signal_members > 1L]),
  n_hidden_variant_gene_pairs = sum(n_hidden_members),
  n_recovered_signal_member_rows = sum(n_signal_members),
  max_members_per_signal = max(n_signal_members),
  median_members_among_tied_signals = if (any(n_signal_members > 1L)) {
    as.numeric(stats::median(n_signal_members[n_signal_members > 1L]))
  } else {
    1.0
  }
), by = .(dataset_id, context, split)]
hidden_unique_by_context <- expanded[was_in_current_pairs == FALSE, .(
  n_hidden_unique_variants = uniqueN(variant_id)
), by = .(dataset_id, context, split)]
summary_by_context <- merge(
  summary_by_context,
  hidden_unique_by_context,
  by = c("dataset_id", "context", "split"),
  all.x = TRUE,
  sort = FALSE
)
summary_by_context[is.na(n_hidden_unique_variants), n_hidden_unique_variants := 0L]
summary_by_context[, pct_row_increase := 100 * (n_recovered_signal_member_rows - n_current_signal_rows) / n_current_signal_rows]

total_summary <- signal_member_counts[, .(
  dataset_id = "ALL",
  context = "ALL",
  split = "all",
  n_current_signal_rows = .N,
  n_signals_with_exact_ties = sum(n_signal_members > 1L),
  n_genes_with_exact_ties = 0L,
  n_hidden_variant_gene_pairs = sum(n_hidden_members),
  n_hidden_unique_variants = uniqueN(expanded[was_in_current_pairs == FALSE, variant_id]),
  n_recovered_signal_member_rows = sum(n_signal_members),
  max_members_per_signal = max(n_signal_members),
  median_members_among_tied_signals = if (any(n_signal_members > 1L)) {
    as.numeric(stats::median(n_signal_members[n_signal_members > 1L]))
  } else {
    1.0
  }
)]
total_summary[, n_genes_with_exact_ties := uniqueN(
  signal_member_counts[n_signal_members > 1L, paste(dataset_id, gene_id, sep = "|")]
)]
total_summary[, pct_row_increase := 100 * (n_recovered_signal_member_rows - n_current_signal_rows) / n_current_signal_rows]
summary_by_context <- rbindlist(list(summary_by_context, total_summary), use.names = TRUE, fill = TRUE)
summary_by_context[, `:=`(
  total_order_tmp = as.integer(context == "ALL"),
  context_order_tmp = match(context, names(SEURAT_CONTEXT_TO_DATASET_ID))
)]
setorder(summary_by_context, total_order_tmp, context_order_tmp)
summary_by_context[, c("total_order_tmp", "context_order_tmp") := NULL]

gwas_summary_parts <- list()
new_gwas_parts <- list()
gwas_specs <- list(
  list(disorder = "MDD", threshold = "strict", flag = "MDD_gwasVar_strict"),
  list(disorder = "MDD", threshold = "exploratory", flag = "MDD_gwasVar_exp"),
  list(disorder = "BD", threshold = "strict", flag = "BD_gwasVar_strict"),
  list(disorder = "BD", threshold = "exploratory", flag = "BD_gwasVar_exp"),
  list(disorder = "SCZD", threshold = "strict", flag = "SCZD_gwasVar_strict")
)

for (spec in gwas_specs) {
  if (!spec$flag %in% names(expanded)) next
  current_signal <- expanded[is_index_variant == TRUE, .(
    current_overlap = any(get(spec$flag) == 1L)
  ), by = .(dataset_id, context, split, gene_id, gene_name, signal_id, index_variant_id)]
  recovered_signal <- expanded[, .(
    recovered_overlap = any(get(spec$flag) == 1L),
    recovered_variants = paste(sort(unique(variant_id[get(spec$flag) == 1L])), collapse = ", ")
  ), by = .(dataset_id, context, split, gene_id, gene_name, signal_id, index_variant_id)]
  cmp <- merge(current_signal, recovered_signal, by = c(
    "dataset_id", "context", "split", "gene_id", "gene_name", "signal_id", "index_variant_id"
  ))
  cmp[, `:=`(
    disorder = spec$disorder,
    threshold = spec$threshold,
    newly_recovered_overlap = current_overlap == FALSE & recovered_overlap == TRUE
  )]
  gwas_summary_parts[[paste(spec$disorder, spec$threshold)]] <- cmp[, .(
    n_signals = .N,
    n_current_exact_overlap_signals = sum(current_overlap),
    n_recovered_exact_overlap_signals = sum(recovered_overlap),
    n_new_exact_overlap_signals = sum(newly_recovered_overlap),
    n_new_exact_overlap_genes = uniqueN(gene_id[newly_recovered_overlap])
  ), by = .(dataset_id, context, split, disorder, threshold)]
  new_gwas_parts[[paste(spec$disorder, spec$threshold)]] <- cmp[newly_recovered_overlap == TRUE]
}

gwas_summary <- rbindlist(gwas_summary_parts, use.names = TRUE, fill = TRUE)
new_gwas_matches <- rbindlist(new_gwas_parts, use.names = TRUE, fill = TRUE)
gwas_summary[, context_order_tmp := match(context, names(SEURAT_CONTEXT_TO_DATASET_ID))]
setorder(gwas_summary, disorder, threshold, context_order_tmp)
gwas_summary[, context_order_tmp := NULL]

coloc <- fread(input_paths$coloc_gated, sep = "\t")
assert_required_cols(
  coloc,
  c("disorder", "dataset_id", "context", "gene_id", "gene_name", "lead_snp", "gate_reason"),
  basename(input_paths$coloc_gated)
)
coloc <- coloc[gate_reason == "pass"]
coloc[, lead_chr := variant_chr(lead_snp)]
coloc <- coloc[dataset_id %in% selected_datasets & lead_chr %in% selected_chromosomes]
coloc[, coloc_row_id := .I]
coloc[, `:=`(current_index_match = FALSE, recovered_signal_match = FALSE)]

current_coloc_lut <- unique(signals[, .(
  dataset_id,
  context,
  gene_id,
  index_variant_id = variant_id
)])
recovered_coloc_lut <- unique(expanded[, .(
  dataset_id,
  context,
  gene_id,
  recovered_variant_id = variant_id,
  index_variant_id,
  signal_id
)])

coloc[current_coloc_lut, on = .(
  dataset_id,
  context,
  gene_id,
  lead_snp = index_variant_id
), current_index_match := TRUE]
coloc[recovered_coloc_lut, on = .(
  dataset_id,
  context,
  gene_id,
  lead_snp = recovered_variant_id
), recovered_signal_match := TRUE]
coloc[, newly_recovered_match := current_index_match == FALSE & recovered_signal_match == TRUE]

coloc_summary <- coloc[, .(
  n_coloc_rows = .N,
  n_current_index_matches = sum(current_index_match),
  n_recovered_signal_matches = sum(recovered_signal_match),
  n_new_signal_matches = sum(newly_recovered_match),
  n_new_match_genes = uniqueN(gene_id[newly_recovered_match])
), by = .(disorder, dataset_id, context)]
new_coloc_matches <- coloc[newly_recovered_match == TRUE]
coloc_summary[, context_order_tmp := match(context, names(SEURAT_CONTEXT_TO_DATASET_ID))]
setorder(coloc_summary, disorder, context_order_tmp)
coloc_summary[, context_order_tmp := NULL]

output_paths <- list(
  recovered_pairs = file.path(output_dir, "map_significant_pairs_tie_recovered.csv.gz"),
  recovered_independent = file.path(output_dir, "map_independent_significant_tie_recovered.csv.gz"),
  members = file.path(output_dir, "independent_eqtl_tie_recovery_members.csv.gz"),
  summary = file.path(output_dir, "independent_eqtl_tie_recovery_summary_by_context.csv"),
  gwas_summary = file.path(output_dir, "independent_eqtl_tie_recovery_gwas_changes_by_context.csv"),
  new_gwas = file.path(output_dir, "independent_eqtl_tie_recovery_new_gwas_matches.csv"),
  coloc_summary = file.path(output_dir, "independent_eqtl_tie_recovery_coloc_changes_by_context.csv"),
  new_coloc = file.path(output_dir, "independent_eqtl_tie_recovery_new_coloc_matches.csv"),
  rejected = file.path(output_dir, "independent_eqtl_tie_recovery_rejected_same_p_candidates.csv")
)

fwrite(expanded, output_paths$recovered_pairs)
fwrite(independent_expanded, output_paths$recovered_independent)
fwrite(
  expanded[, intersect(c(metadata_front, "DEG", "MDD_DEG", "BD_DEG", grep("_gwas", names(expanded), value = TRUE),
                          "member_raw_pval_nominal", "member_raw_slope", "member_raw_slope_se"), names(expanded)), with = FALSE],
  output_paths$members
)
fwrite(summary_by_context, output_paths$summary)
fwrite(gwas_summary, output_paths$gwas_summary)
fwrite(new_gwas_matches, output_paths$new_gwas)
fwrite(coloc_summary, output_paths$coloc_summary)
fwrite(new_coloc_matches, output_paths$new_coloc)
fwrite(rejected_same_p, output_paths$rejected)

if (isTRUE(opts$promote_current_pairs)) {
  index_backup <- sub("[.]csv[.]gz$", "_index_only.csv.gz", input_paths$current_pairs)
  if (identical(index_backup, input_paths$current_pairs)) {
    index_backup <- paste0(input_paths$current_pairs, ".index_only.csv.gz")
  }
  fwrite(current_pairs, index_backup)
  fwrite(expanded, input_paths$current_pairs)
  message("Promoted recovered pairs to canonical current-pairs table: ", input_paths$current_pairs)
  message("Preserved the one-index recovery baseline: ", index_backup)
}

message("Recovery complete. Outputs:")
for (path in output_paths) message("  ", path)
