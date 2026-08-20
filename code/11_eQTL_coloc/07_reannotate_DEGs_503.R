#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})

args <- commandArgs(trailingOnly = TRUE)
install_outputs <- "--install" %in% args
args <- setdiff(args, "--install")

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Cannot determine script path")
script_file <- normalizePath(sub("^--file=", "", script_arg))
code_dir <- dirname(script_file)
repo_root <- normalizePath(file.path(code_dir, "..", ".."))
source(file.path(code_dir, "utils.R"), chdir = FALSE)

default_output_dir <- file.path(code_dir, "outputs", "jacqui_deg503_20260819", "candidate")
output_dir <- if (length(args) == 0L) default_output_dir else normalizePath(args[[1]], mustWork = FALSE)
tables_dir <- file.path(repo_root, "processed-data", "11_eQTL_coloc", "seurat", "tables")
final_dir <- file.path(repo_root, "processed-data", "11_eQTL_coloc", "seurat", "final")
snapshot_dir <- file.path(code_dir, "outputs", "jacqui_deg503_20260819", "pre_change")
source_tables_dir <- file.path(snapshot_dir, "tables")
source_final_dir <- file.path(snapshot_dir, "final")
candidate_tables_dir <- file.path(output_dir, "tables")
candidate_final_dir <- file.path(output_dir, "final")
report_dir <- file.path(output_dir, "reports")

for (path in c(candidate_tables_dir, candidate_final_dir, report_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

fail <- function(...) stop(..., call. = FALSE)

assert_true <- function(value, message) {
  if (!isTRUE(value)) fail(message)
  invisible(TRUE)
}

assert_columns <- function(dt, columns, label) {
  missing <- setdiff(columns, names(dt))
  if (length(missing) > 0L) {
    fail(label, " is missing columns: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

same_values <- function(left, right) {
  isTRUE(all.equal(left, right, check.attributes = FALSE))
}

same_text <- function(left, right) {
  left <- as.character(left)
  right <- as.character(right)
  left[is.na(left)] <- ""
  right[is.na(right)] <- ""
  identical(left, right)
}

same_serialized_values <- function(left, right) {
  if (!identical(names(left), names(right)) || nrow(left) != nrow(right)) return(FALSE)
  checks <- vapply(names(left), function(column) {
    if (is.character(left[[column]]) || is.character(right[[column]])) {
      same_text(left[[column]], right[[column]])
    } else {
      same_values(left[[column]], right[[column]])
    }
  }, logical(1))
  all(checks)
}

read_table <- function(path) {
  if (!file.exists(path)) fail("Missing required table: ", path)
  fread(path)
}

write_checked <- function(dt, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(dt, path)
  observed <- fread(path)
  if (!same_serialized_values(dt, observed)) fail("Round-trip validation failed: ", path)
  invisible(path)
}

write_deg_only_preserving_text <- function(source_path, dt, path) {
  raw <- fread(source_path, colClasses = "character", na.strings = NULL)
  if (!identical(names(raw), names(dt)) || nrow(raw) != nrow(dt)) {
    fail("Text-preserving source does not match typed table: ", source_path)
  }
  non_deg_cols <- setdiff(names(raw), "DEG")
  raw[, DEG := as.character(dt$DEG)]
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(raw, path)
  observed_raw <- fread(path, colClasses = "character", na.strings = NULL)
  if (!identical(raw[, ..non_deg_cols], observed_raw[, ..non_deg_cols])) {
    fail("Non-DEG text changed while writing ", path)
  }
  if (!identical(raw$DEG, observed_raw$DEG)) fail("DEG text changed while writing ", path)
  observed <- fread(path)
  if (!same_serialized_values(dt, observed)) fail("Typed round-trip validation failed: ", path)
  invisible(path)
}

atomic_replace <- function(source, target) {
  if (!file.exists(source)) fail("Missing validated candidate: ", source)
  temp_target <- paste0(target, ".deg503.tmp")
  if (file.exists(temp_target)) unlink(temp_target)
  if (!file.copy(source, temp_target, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)) {
    fail("Failed to prepare atomic replacement for ", target)
  }
  status <- system2("mv", c("-f", temp_target, target))
  if (!identical(status, 0L)) fail("Failed to atomically replace ", target)
  invisible(target)
}

snapshot_targets <- c(
  setNames(
    file.path(tables_dir, list.files(file.path(snapshot_dir, "tables"))),
    file.path(snapshot_dir, "tables", list.files(file.path(snapshot_dir, "tables")))
  ),
  setNames(
    file.path(final_dir, list.files(file.path(snapshot_dir, "final"))),
    file.path(snapshot_dir, "final", list.files(file.path(snapshot_dir, "final")))
  )
)
workspace_snapshot_matches <- 0L
workspace_candidate_matches <- 0L
for (snapshot_file in names(snapshot_targets)) {
  current_file <- snapshot_targets[[snapshot_file]]
  if (!file.exists(current_file)) fail("Current output missing for snapshot file: ", current_file)
  candidate_group <- basename(dirname(snapshot_file))
  candidate_file <- file.path(output_dir, candidate_group, basename(snapshot_file))
  comparison_files <- c(snapshot_file, current_file)
  if (file.exists(candidate_file)) comparison_files <- c(comparison_files, candidate_file)
  hashes <- unname(tools::md5sum(comparison_files))
  matches_snapshot <- identical(hashes[[1]], hashes[[2]])
  matches_candidate <- length(hashes) == 3L && identical(hashes[[2]], hashes[[3]])
  if (!matches_snapshot && !matches_candidate) {
    fail("Workspace output contradicts the pre-change snapshot: ", current_file)
  }
  if (matches_snapshot) {
    workspace_snapshot_matches <- workspace_snapshot_matches + 1L
  } else {
    workspace_candidate_matches <- workspace_candidate_matches + 1L
  }
}

degs <- load_DEGs(repo_root = repo_root, mode = "standard", verbose = FALSE)
legacy_degs <- load_legacy_standard_DEGs(repo_root = repo_root, verbose = FALSE)
sig_df <- as.data.table(degs$tables$sig_df)
deg_global <- as.data.table(degs$global)
authoritative_ids <- unique(deg_global$gene_id)
authoritative_names <- unique(deg_global$gene_name)
legacy_ids <- unique(as.data.table(legacy_degs$global)$gene_id)
disorder_ids <- lapply(degs$disorder_related$global, function(x) unique(as.data.table(x)$gene_id))

assert_true(nrow(sig_df) == 2140L, "Jacqui sig.df does not contain 2,140 rows")
assert_true(length(authoritative_ids) == 503L, "Authoritative DEG set does not contain 503 gene IDs")
assert_true(length(authoritative_names) == 503L, "Authoritative DEG set does not contain 503 gene names")
assert_true(length(legacy_ids) == 817L, "Legacy DEG set does not contain 817 gene IDs")
assert_true(length(setdiff(legacy_ids, authoritative_ids)) == 314L, "Legacy-to-authoritative removal is not 314 genes")
assert_true(length(setdiff(authoritative_ids, legacy_ids)) == 0L, "Authoritative DEG set adds genes outside the legacy set")
assert_true(length(disorder_ids$MDD) == 224L, "Authoritative MDD DEG set does not contain 224 genes")
assert_true(length(disorder_ids$BD) == 281L, "Authoritative BD DEG set does not contain 281 genes")

row_files <- c(
  "map_cis_significant.csv.gz",
  "map_independent_significant.csv.gz",
  "map_significant_unified.csv.gz",
  "map_significant_pairs.csv.gz",
  "map_significant_pairs_index_only.csv.gz",
  "map_significant_pairs_tie_recovered.csv.gz",
  "map_independent_significant_tie_recovered.csv.gz",
  "independent_eqtl_tie_recovery_members.csv.gz",
  "nominal_BH05.csv.gz"
)
summary_files <- c(
  "map_cis_summary.csv",
  "map_independent_summary.csv",
  "map_significant_summary.csv",
  "nominal_BH05_summary.csv"
)
helper_files <- c("eqtl_boxplot_deg_pairs.csv", "eqtl_boxplot_examples.csv")
slide_file <- "slide_table_seurat_sczd_gwas.csv"
handled_files <- c(row_files, summary_files, helper_files, slide_file)

table_candidates <- list.files(tables_dir, pattern = "[.]csv([.]gz)?$", full.names = TRUE)
for (path in table_candidates) {
  header <- names(fread(path, nrows = 0L))
  relevant <- grepl("(^DEG$|^n_DEG$|^DEG_genes$|trifecta|DEGxGWAS|^DEGs$|^nDEG)", header)
  if (any(relevant) && !basename(path) %in% handled_files) {
    fail("Unrecognized table with DEG-derived columns: ", path)
  }
}

projected <- data.table(
  table = c(
    "map_cis_significant.csv.gz",
    "map_independent_significant.csv.gz",
    "map_significant_unified.csv.gz",
    "map_significant_pairs.csv.gz",
    "map_significant_pairs_tie_recovered.csv.gz",
    "map_independent_significant_tie_recovered.csv.gz",
    "nominal_BH05.csv.gz"
  ),
  rows = c(3296L, 3405L, 6701L, 15677L, 15677L, 15591L, 455143L),
  old_deg_rows = c(242L, 258L, 500L, 1133L, 1133L, 1123L, 40766L),
  new_deg_rows = c(138L, 143L, 281L, 473L, 473L, 465L, 10420L),
  removed_rows = c(104L, 115L, 219L, 660L, 660L, 658L, 30346L),
  old_deg_genes = c(77L, 77L, 77L, 77L, 77L, 77L, 316L),
  new_deg_genes = c(50L, 50L, 50L, 50L, 50L, 50L, 188L)
)

old_tables <- list()
new_tables <- list()
impact_tables <- list()
impact_contexts <- list()
impact_overlaps <- list()

for (file_name in row_files) {
  source_path <- file.path(source_tables_dir, file_name)
  old <- read_table(source_path)
  assert_columns(old, c("gene_id", "gene_name", "DEG"), file_name)
  new <- copy(old)
  new[, DEG := as.integer(gene_id %in% authoritative_ids)]

  if (any(old$DEG == 0L & new$DEG == 1L, na.rm = TRUE)) {
    fail("Unexpected 0 -> 1 DEG change in ", file_name)
  }
  non_deg_cols <- setdiff(names(old), "DEG")
  if (!same_values(old[, ..non_deg_cols], new[, ..non_deg_cols])) {
    fail("Non-DEG values changed in ", file_name)
  }
  for (disorder in intersect(c("MDD", "BD"), sub("_DEG$", "", grep("_(DEG)$", names(old), value = TRUE)))) {
    flag_col <- paste0(disorder, "_DEG")
    if (!identical(old[[flag_col]], new[[flag_col]])) fail(flag_col, " changed in ", file_name)
    expected_flag <- as.integer(old$gene_id %in% disorder_ids[[disorder]])
    if (!identical(as.integer(old[[flag_col]]), expected_flag)) {
      fail(flag_col, " does not match authoritative sig.df in ", file_name)
    }
  }

  table_impact <- data.table(
    table = file_name,
    rows = nrow(old),
    old_deg_rows = sum(old$DEG == 1L, na.rm = TRUE),
    new_deg_rows = sum(new$DEG == 1L, na.rm = TRUE),
    removed_rows = sum(old$DEG == 1L & new$DEG == 0L, na.rm = TRUE),
    added_rows = sum(old$DEG == 0L & new$DEG == 1L, na.rm = TRUE),
    old_deg_genes = uniqueN(old[DEG == 1L, gene_id]),
    new_deg_genes = uniqueN(new[DEG == 1L, gene_id])
  )
  impact_tables[[file_name]] <- table_impact

  if ("context" %in% names(old)) {
    impact_contexts[[file_name]] <- old[, .(
      rows = .N,
      old_deg_rows = sum(DEG == 1L, na.rm = TRUE),
      old_deg_genes = uniqueN(gene_id[DEG == 1L])
    ), by = context][new[, .(
      context,
      new_deg_rows = sum(DEG == 1L, na.rm = TRUE),
      new_deg_genes = uniqueN(gene_id[DEG == 1L])
    ), by = context], on = "context"][, `:=`(
      table = file_name,
      removed_rows = old_deg_rows - new_deg_rows,
      removed_genes = old_deg_genes - new_deg_genes
    )]
  }

  for (disorder in c("MDD", "BD", "SCZD")) {
    for (level in c("strict", "exp")) {
      flag_col <- paste0(disorder, "_gwas_", level)
      if (!flag_col %in% names(old)) next
      impact_overlaps[[paste(file_name, disorder, level, sep = "|")]] <- data.table(
        table = file_name,
        disorder = disorder,
        level = level,
        old_deg_gwas_rows = sum(old$DEG == 1L & old[[flag_col]] == 1L, na.rm = TRUE),
        new_deg_gwas_rows = sum(new$DEG == 1L & new[[flag_col]] == 1L, na.rm = TRUE),
        old_deg_gwas_genes = uniqueN(old[DEG == 1L & get(flag_col) == 1L, gene_id]),
        new_deg_gwas_genes = uniqueN(new[DEG == 1L & get(flag_col) == 1L, gene_id])
      )
    }
  }

  expected <- projected[table == file_name]
  if (nrow(expected) == 1L) {
    observed <- table_impact[, .(rows, old_deg_rows, new_deg_rows, removed_rows, old_deg_genes, new_deg_genes)]
    target <- expected[, .(rows, old_deg_rows, new_deg_rows, removed_rows, old_deg_genes, new_deg_genes)]
    if (!identical(observed, target)) {
      print(observed)
      print(target)
      fail("Projected annotation impact mismatch for ", file_name)
    }
  }

  write_deg_only_preserving_text(source_path, new, file.path(candidate_tables_dir, file_name))
  old_tables[[file_name]] <- old
  new_tables[[file_name]] <- new
}

canonical <- new_tables[["map_significant_pairs.csv.gz"]]
tie_recovered <- new_tables[["map_significant_pairs_tie_recovered.csv.gz"]]
pair_key <- c("dataset_id", "context", "split", "gene_id", "variant_id")
assert_true(nrow(canonical) == 15677L, "Canonical significant-pair table does not contain 15,677 rows")
assert_true(uniqueN(canonical$signal_id) == 3435L, "Canonical significant-pair table does not contain 3,435 signals")
assert_true(canonical[is_index_variant == TRUE, .N] == 3435L, "Canonical table does not contain 3,435 representative rows")
assert_true(canonical[recovery_status == "recovered_exact_tie", .N] == 12242L, "Canonical table does not contain 12,242 recovered ties")
assert_true(anyDuplicated(canonical[, ..pair_key]) == 0L, "Canonical pair keys are duplicated")
assert_true(same_values(canonical, tie_recovered), "Canonical and tie-recovered pair tables differ after reannotation")

summary_sources <- list(
  "map_cis_summary.csv" = canonical[cis_supported == TRUE],
  "map_independent_summary.csv" = canonical[indep_supported == TRUE],
  "map_significant_summary.csv" = canonical,
  "nominal_BH05_summary.csv" = new_tables[["nominal_BH05.csv.gz"]]
)
new_summaries <- list()
summary_impacts <- list()

for (file_name in names(summary_sources)) {
  old <- read_table(file.path(source_tables_dir, file_name))
  source_dt <- summary_sources[[file_name]]
  assert_columns(old, c("split", "context", "n_eGenes", "n_DEG", "DEG_genes"), file_name)
  new <- copy(old)

  for (row_index in seq_len(nrow(new))) {
    split_value <- new$split[[row_index]]
    context_value <- new$context[[row_index]]
    group <- source_dt[split == split_value & context == context_value]
    if (uniqueN(group$gene_id) != old$n_eGenes[[row_index]]) {
      fail("Unchanged eGene count does not reconcile for ", file_name, " / ", split_value, " / ", context_value)
    }
    new$n_DEG[[row_index]] <- uniqueN(group[DEG == 1L, gene_id])
    new$DEG_genes[[row_index]] <- collapse_gene_list(group[DEG == 1L, gene_name])

    trifecta_cols <- grep("^n_trifecta_", names(new), value = TRUE)
    for (count_col in trifecta_cols) {
      tokens <- strsplit(sub("^n_trifecta_", "", count_col), "_", fixed = TRUE)[[1]]
      disorder <- tokens[[1]]
      level <- tokens[[2]]
      flag_col <- paste0(disorder, "_gwas_", level)
      if (!flag_col %in% names(group)) fail("Missing ", flag_col, " for ", file_name)
      keep <- group$DEG == 1L & group[[flag_col]] == 1L
      new[[count_col]][[row_index]] <- uniqueN(group$gene_id[keep])
      gene_col <- sub("^n_", "", paste0(count_col, "_genes"))
      if (gene_col %in% names(new)) new[[gene_col]][[row_index]] <- collapse_gene_list(group$gene_name[keep])
    }

    for (disorder in c("MDD", "BD")) {
      disorder_flag <- paste0(disorder, "_DEG")
      disorder_count <- paste0("n_", disorder_flag)
      disorder_genes <- paste0(disorder_flag, "_genes")
      if (disorder_count %in% names(new)) {
        expected_flag <- as.integer(group$gene_id %in% disorder_ids[[disorder]])
        if (disorder_flag %in% names(group) && !identical(as.integer(group[[disorder_flag]]), expected_flag)) {
          fail(disorder_flag, " changed or is inconsistent in source for ", file_name)
        }
        expected_count <- uniqueN(group$gene_id[expected_flag == 1L])
        expected_genes <- collapse_gene_list(group$gene_name[expected_flag == 1L])
        if (old[[disorder_count]][[row_index]] != expected_count) fail(disorder_count, " unexpectedly differs in ", file_name)
        if (!same_text(old[[disorder_genes]][[row_index]], expected_genes)) fail(disorder_genes, " unexpectedly differs in ", file_name)
      }

      for (level in c("strict", "exp")) {
        count_col <- paste0("n_", disorder, "_DEGxGWAS_", level)
        gene_col <- paste0(disorder, "_DEGxGWAS_", level, "_genes")
        flag_col <- paste0(disorder, "_gwas_", level)
        if (!count_col %in% names(new)) next
        expected_flag <- as.integer(group$gene_id %in% disorder_ids[[disorder]])
        keep <- expected_flag == 1L & group[[flag_col]] == 1L
        if (old[[count_col]][[row_index]] != uniqueN(group$gene_id[keep])) fail(count_col, " unexpectedly differs in ", file_name)
        if (!same_text(old[[gene_col]][[row_index]], collapse_gene_list(group$gene_name[keep]))) fail(gene_col, " unexpectedly differs in ", file_name)
      }
    }
  }

  general_derived <- c(
    "n_DEG", "DEG_genes",
    grep("^(n_trifecta_|trifecta_.*_genes$)", names(new), value = TRUE)
  )
  unchanged_cols <- setdiff(names(old), general_derived)
  if (!same_values(old[, ..unchanged_cols], new[, ..unchanged_cols])) {
    fail("Non-general-DEG summary values changed in ", file_name)
  }
  summary_impacts[[file_name]] <- new[, .(
    table = file_name,
    split,
    context,
    old_n_DEG = old$n_DEG,
    new_n_DEG = n_DEG,
    removed_DEG_eGenes = old$n_DEG - n_DEG
  )]
  write_checked(new, file.path(candidate_tables_dir, file_name))
  new_summaries[[file_name]] <- new
}

boxplot_old <- read_table(file.path(source_tables_dir, "eqtl_boxplot_deg_pairs.csv"))
assert_columns(boxplot_old, c("gene_id", "gene_name", "DEG"), "eqtl_boxplot_deg_pairs.csv")
boxplot_new <- copy(boxplot_old[gene_id %in% authoritative_ids])
boxplot_new[, DEG := 1L]
assert_true(all(boxplot_new$DEG == 1L), "Filtered boxplot DEG metadata contains DEG == 0")
assert_true(anyDuplicated(boxplot_new[, .(dataset_id, gene_id, variant_id)]) == 0L, "Filtered boxplot DEG metadata has duplicate keys")
write_checked(boxplot_new, file.path(candidate_tables_dir, "eqtl_boxplot_deg_pairs.csv"))

examples_old <- read_table(file.path(source_tables_dir, "eqtl_boxplot_examples.csv"))
assert_columns(examples_old, c("gene_id", "DEG"), "eqtl_boxplot_examples.csv")
examples_new <- copy(examples_old)
examples_new[, DEG := as.integer(gene_id %in% authoritative_ids)]
assert_true(!any(examples_old$DEG == 0L & examples_new$DEG == 1L), "Unexpected 0 -> 1 change in boxplot examples")
write_checked(examples_new, file.path(candidate_tables_dir, "eqtl_boxplot_examples.csv"))

slide_old <- read_table(file.path(source_tables_dir, slide_file))
assert_columns(slide_old, c("DEGs", "nDEG_nGWAS", "nDEGS", "nDEG_nGWAS_genes"), slide_file)
split_gene_list <- function(value) {
  if (is.na(value) || !nzchar(value)) return(character())
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}
collapse_names <- function(values) paste(sort(unique(values[nzchar(values)])), collapse = ", ")
slide_new <- copy(slide_old)
for (row_index in seq_len(nrow(slide_new))) {
  deg_names <- intersect(split_gene_list(slide_old$nDEGS[[row_index]]), authoritative_names)
  overlap_names <- intersect(split_gene_list(slide_old$nDEG_nGWAS_genes[[row_index]]), authoritative_names)
  slide_new$DEGs[[row_index]] <- length(unique(deg_names))
  slide_new$nDEGS[[row_index]] <- collapse_names(deg_names)
  slide_new$nDEG_nGWAS[[row_index]] <- length(unique(overlap_names))
  slide_new$nDEG_nGWAS_genes[[row_index]] <- collapse_names(overlap_names)
}
slide_unchanged <- setdiff(names(slide_old), c("DEGs", "nDEG_nGWAS", "nDEGS", "nDEG_nGWAS_genes"))
assert_true(same_values(slide_old[, ..slide_unchanged], slide_new[, ..slide_unchanged]), "Non-DEG slide-table values changed")
write_checked(slide_new, file.path(candidate_tables_dir, slide_file))

context_order <- names(SEURAT_CONTEXT_TO_DATASET_ID)
dataset_order <- unname(SEURAT_CONTEXT_TO_DATASET_ID)
split_order <- DEG_SPLITS

set_context_order <- function(dt, extra_order = character()) {
  rank_cols <- character()
  if ("context" %in% names(dt)) {
    dt[, context_rank_tmp := match(context, context_order)]
    rank_cols <- c(rank_cols, "context_rank_tmp")
  }
  if ("dataset_id" %in% names(dt)) {
    dt[, dataset_rank_tmp := match(dataset_id, dataset_order)]
    rank_cols <- c(rank_cols, "dataset_rank_tmp")
  }
  if ("split" %in% names(dt)) {
    dt[, split_rank_tmp := match(split, split_order)]
    rank_cols <- c(rank_cols, "split_rank_tmp")
  }
  order_cols <- intersect(c(extra_order, rank_cols, "gene_name", "gene_id", "variant_id"), names(dt))
  if (length(order_cols) > 0L) setorderv(dt, order_cols)
  if (length(rank_cols) > 0L) dt[, (rank_cols) := NULL]
  invisible(dt)
}

front_order <- function(dt, front_cols) {
  setcolorder(dt, c(intersect(front_cols, names(dt)), setdiff(names(dt), front_cols)))
  invisible(dt)
}

variant_info <- load_genotype_variant_info(
  plink2_prefix = file.path(repo_root, "processed-data", "00_genotypes", "plink2", "merged_maf05"),
  use_cache = TRUE,
  build_if_missing = FALSE
)
rsid_lut <- variant_info[, .(variant_id, rsid)]
setkey(rsid_lut, variant_id)

add_rsids <- function(dt, label) {
  out <- copy(dt)
  out[rsid_lut, on = "variant_id", rsid := i.rsid]
  if (anyNA(out$rsid)) fail(label, " has variants missing from the existing genotype variant cache")
  out
}

cis_public_source <- add_rsids(canonical, "cis_independent")
nominal_public_source <- add_rsids(new_tables[["nominal_BH05.csv.gz"]], "nominal_BH05")

make_cis_public <- function(dt) {
  out <- copy(dt)
  out[, eQTL_signal_class := fcase(
    pair_provenance == "cis_and_independent", "lead_cis_and_independent",
    pair_provenance == "cis_only", "lead_cis_only",
    pair_provenance == "independent_only", "independent_signal",
    default = pair_provenance
  )]
  out[, lead_cis_qval := fifelse(!is.na(qval), qval, qval_parent)]
  out[, independent_signal_rank := fifelse(
    signal_kind == "conditionally_independent",
    as.integer(signal_rank),
    NA_integer_
  )]
  out[, eQTL_signal_id := signal_id]
  out[!nzchar(rsid), rsid := NA_character_]
  columns <- c(
    "context", "gene_id", "gene_name", "eQTL_signal_id", "eQTL_signal_class",
    "independent_signal_rank", "variant_id", "rsid", "DEG", "MDD_DEG", "BD_DEG",
    "MDD_gwasVar_strict", "MDD_gwasVar_exp", "MDD_gwasGene", "MDD_gwas_strict", "MDD_gwas_exp",
    "BD_gwasVar_strict", "BD_gwasVar_exp", "BD_gwasGene", "BD_gwas_strict", "BD_gwas_exp",
    "SCZD_gwasVar_strict", "SCZD_gwasGene", "SCZD_gwas_strict",
    "start_distance", "end_distance", "af", "pval_nominal", "slope", "slope_se", "pval_perm", "lead_cis_qval"
  )
  assert_columns(out, columns, "public cis_independent")
  out <- out[, ..columns]
  set_context_order(out, extra_order = c("context_rank_tmp", "gene_name", "pval_nominal"))
  out
}

make_nominal_public <- function(dt) {
  out <- copy(dt)
  out[!nzchar(rsid), rsid := NA_character_]
  columns <- c(
    "context", "gene_id", "gene_name", "variant_id", "rsid", "DEG",
    "MDD_gwasVar_strict", "MDD_gwasGene", "MDD_gwas_strict",
    "BD_gwasVar_strict", "BD_gwasGene", "BD_gwas_strict",
    "SCZD_gwasVar_strict", "SCZD_gwasGene", "SCZD_gwas_strict",
    "start_distance", "af", "pval_nominal", "fdr", "slope", "slope_se"
  )
  assert_columns(out, columns, "public nominal_BH05")
  out <- out[, ..columns]
  set_context_order(out, extra_order = c("context_rank_tmp", "gene_name", "fdr", "pval_nominal"))
  out
}

summary_source <- function(dt, eqtl_type_value) {
  out <- copy(dt)
  out[, dataset_id := unname(SEURAT_CONTEXT_TO_DATASET_ID[context])]
  out[, eQTL_type := eqtl_type_value]
  front_order(out, c("eQTL_type", "dataset_id", "context", "split"))
  out
}

manifest_qc <- read_table(file.path(source_tables_dir, "manifest_qc.csv"))
if (!"context" %in% names(manifest_qc) && "seurat_label" %in% names(manifest_qc)) {
  manifest_qc[, context := seurat_label]
}
eqtl_summary <- rbindlist(
  list(
    summary_source(new_summaries[["map_cis_summary.csv"]], "lead_cis_qval05"),
    summary_source(new_summaries[["map_independent_summary.csv"]], "independent_parent_qval05_perm05"),
    summary_source(new_summaries[["map_significant_summary.csv"]], "cis_independent"),
    summary_source(new_summaries[["nominal_BH05_summary.csv"]], "nominal_BH05")
  ),
  use.names = TRUE,
  fill = TRUE
)
manifest_context <- unique(manifest_qc[, .(
  dataset_id, context, split, covariate_model, n_samples, n_genes_bed, n_genes_pca, n_expr_pcs
)])
eqtl_summary <- merge(
  eqtl_summary,
  manifest_context,
  by = c("dataset_id", "context", "split"),
  all.x = TRUE,
  sort = FALSE
)
front_order(eqtl_summary, c(
  "eQTL_type", "dataset_id", "context", "split", "covariate_model",
  "n_samples", "n_genes_bed", "n_genes_pca", "n_expr_pcs"
))
eqtl_summary[, eQTL_type_rank_tmp := match(
  eQTL_type,
  c("lead_cis_qval05", "independent_parent_qval05_perm05", "cis_independent", "nominal_BH05")
)]
set_context_order(eqtl_summary, extra_order = c("eQTL_type_rank_tmp"))
eqtl_summary[, eQTL_type_rank_tmp := NULL]

make_eqtl_summary_public <- function(dt) {
  out <- copy(dt)
  if ("n_genes_bed" %in% names(out)) setnames(out, "n_genes_bed", "n_genes_tested")
  count_cols <- grep("^n_", names(out), value = TRUE)
  drop_counts <- c("n_genes_pca", "n_expr_pcs")
  columns <- c(
    "eQTL_type", "context", "n_samples", "n_genes_tested",
    setdiff(count_cols, c("n_samples", "n_genes_tested", drop_counts))
  )
  out[, ..columns]
}

max_flag <- function(values) {
  values <- as.integer(values)
  if (all(is.na(values))) return(NA_integer_)
  as.integer(max(values, na.rm = TRUE))
}

collapse_gene_list_limited <- function(values, max_chars = 30000L) {
  values <- sort(unique(as.character(values[!is.na(values) & nzchar(as.character(values))])))
  if (length(values) == 0L) return(list(text = "", n_total = 0L, n_listed = 0L, truncated = FALSE))
  text <- ""
  n_listed <- 0L
  for (value in values) {
    candidate <- if (nzchar(text)) paste(text, value, sep = ", ") else value
    if (nchar(candidate, type = "chars", allowNA = FALSE) > max_chars) break
    text <- candidate
    n_listed <- n_listed + 1L
  }
  truncated <- n_listed < length(values)
  if (truncated) {
    suffix <- sprintf(" ... [truncated; %d of %d listed]", n_listed, length(values))
    available <- max_chars - nchar(suffix, type = "chars", allowNA = FALSE)
    text <- if (available > 0L) paste0(substr(text, 1L, available), suffix) else substr(text, 1L, max_chars)
  }
  list(text = text, n_total = length(values), n_listed = n_listed, truncated = truncated)
}

make_gene_level <- function(dt) {
  flag_cols <- grep("_(DEG|gwasVar_strict|gwasVar_exp|gwasGene|gwas_strict|gwas_exp)$|^DEG$", names(dt), value = TRUE)
  out <- unique(dt[, c("dataset_id", "context", "split", "gene_id", "gene_name", flag_cols), with = FALSE])
  out <- out[, lapply(.SD, max_flag), by = .(dataset_id, context, split, gene_id, gene_name), .SDcols = flag_cols]
  out[, DEG := as.integer(gene_id %in% authoritative_ids)]
  for (disorder in c("MDD", "BD")) out[, (paste0(disorder, "_DEG")) := as.integer(gene_id %in% disorder_ids[[disorder]])]
  out
}

append_gene_row <- function(rows, eqtl_type, subset_dt, disorder, annotation_set, selected) {
  gene_ids <- collapse_gene_list_limited(selected$gene_id)
  gene_names <- collapse_gene_list_limited(selected$gene_name)
  rows[[length(rows) + 1L]] <- data.table(
    eQTL_type = eqtl_type,
    dataset_id = subset_dt$dataset_id[[1]],
    context = subset_dt$context[[1]],
    split = subset_dt$split[[1]],
    disorder = disorder,
    annotation_set = annotation_set,
    n_genes = gene_ids$n_total,
    gene_ids = gene_ids$text,
    gene_ids_listed = gene_ids$n_listed,
    gene_ids_truncated = gene_ids$truncated,
    gene_names = gene_names$text,
    gene_names_listed = gene_names$n_listed,
    gene_names_truncated = gene_names$truncated
  )
  rows
}

build_intersections <- function(dt, eqtl_type) {
  gene_dt <- make_gene_level(dt)
  groups <- unique(gene_dt[, .(dataset_id, context, split)])
  set_context_order(groups)
  rows <- list()
  for (group_index in seq_len(nrow(groups))) {
    key <- groups[group_index]
    subset_dt <- gene_dt[key, on = .(dataset_id, context, split)]
    rows <- append_gene_row(rows, eqtl_type, subset_dt, "all", "all_eGenes", subset_dt)
    rows <- append_gene_row(rows, eqtl_type, subset_dt, "all", "broad_DEG_eGenes", subset_dt[DEG == 1L])
    for (disorder in c("MDD", "BD")) {
      disorder_col <- paste0(disorder, "_DEG")
      rows <- append_gene_row(rows, eqtl_type, subset_dt, disorder, paste0(disorder, "_DEG_eGenes"), subset_dt[get(disorder_col) == 1L])
    }
    for (disorder in c("MDD", "BD", "SCZD")) {
      strict_var <- paste0(disorder, "_gwasVar_strict")
      strict_gene <- paste0(disorder, "_gwasGene")
      strict_any <- paste0(disorder, "_gwas_strict")
      exp_var <- paste0(disorder, "_gwasVar_exp")
      exp_any <- paste0(disorder, "_gwas_exp")
      disorder_deg <- paste0(disorder, "_DEG")
      if (strict_var %in% names(subset_dt)) rows <- append_gene_row(rows, eqtl_type, subset_dt, disorder, paste0(disorder, "_gwasVar_strict_eGenes"), subset_dt[get(strict_var) == 1L])
      if (strict_gene %in% names(subset_dt)) rows <- append_gene_row(rows, eqtl_type, subset_dt, disorder, paste0(disorder, "_gwasGene_eGenes"), subset_dt[get(strict_gene) == 1L])
      if (strict_any %in% names(subset_dt)) {
        rows <- append_gene_row(rows, eqtl_type, subset_dt, disorder, paste0(disorder, "_gwas_strict_eGenes"), subset_dt[get(strict_any) == 1L])
        rows <- append_gene_row(rows, eqtl_type, subset_dt, disorder, paste0("broad_DEG_x_", disorder, "_gwas_strict"), subset_dt[DEG == 1L & get(strict_any) == 1L])
        if (disorder_deg %in% names(subset_dt)) rows <- append_gene_row(rows, eqtl_type, subset_dt, disorder, paste0(disorder, "_DEG_x_", disorder, "_gwas_strict"), subset_dt[get(disorder_deg) == 1L & get(strict_any) == 1L])
      }
      if (exp_var %in% names(subset_dt)) rows <- append_gene_row(rows, eqtl_type, subset_dt, disorder, paste0(disorder, "_gwasVar_exp_eGenes"), subset_dt[get(exp_var) == 1L])
      if (exp_any %in% names(subset_dt)) {
        rows <- append_gene_row(rows, eqtl_type, subset_dt, disorder, paste0(disorder, "_gwas_exp_eGenes"), subset_dt[get(exp_any) == 1L])
        rows <- append_gene_row(rows, eqtl_type, subset_dt, disorder, paste0("broad_DEG_x_", disorder, "_gwas_exp"), subset_dt[DEG == 1L & get(exp_any) == 1L])
        if (disorder_deg %in% names(subset_dt)) rows <- append_gene_row(rows, eqtl_type, subset_dt, disorder, paste0(disorder, "_DEG_x_", disorder, "_gwas_exp"), subset_dt[get(disorder_deg) == 1L & get(exp_any) == 1L])
      }
    }
  }
  out <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  set_context_order(out, extra_order = c("eQTL_type", "disorder", "annotation_set"))
  out
}

cis_public <- make_cis_public(cis_public_source)
nominal_public <- make_nominal_public(nominal_public_source)
eqtl_summary_public <- make_eqtl_summary_public(eqtl_summary)
intersections <- rbindlist(
  list(
    build_intersections(cis_public_source, "cis_independent"),
    build_intersections(nominal_public_source, "nominal_BH05")
  ),
  use.names = TRUE,
  fill = TRUE
)
intersections_public <- intersections[n_genes > 0L, .(
  eQTL_type, context, disorder, annotation_set, n_genes,
  gene_ids, gene_ids_listed, gene_ids_truncated, gene_names
)]
set_context_order(intersections_public, extra_order = c("eQTL_type", "context_rank_tmp", "disorder", "annotation_set"))

coloc_workbook <- file.path(source_final_dir, "coloc_results.xlsx")
if (!requireNamespace("openxlsx", quietly = TRUE)) fail("openxlsx is required to read the existing coloc workbook")
coloc_pass <- as.data.table(openxlsx::read.xlsx(coloc_workbook, sheet = "coloc_pass", na.strings = ""))
coloc_summary <- as.data.table(openxlsx::read.xlsx(coloc_workbook, sheet = "coloc_summary", na.strings = ""))
assert_columns(coloc_pass, c("gene_id", "DEG", "MDD_DEG", "BD_DEG"), "coloc_pass")
old_coloc_deg <- copy(coloc_pass$DEG)
coloc_pass[, DEG := as.integer(gene_id %in% authoritative_ids)]
assert_true(!any(old_coloc_deg == 0L & coloc_pass$DEG == 1L), "Unexpected 0 -> 1 DEG change in coloc_pass")
assert_true(identical(as.integer(coloc_pass$MDD_DEG), as.integer(coloc_pass$gene_id %in% disorder_ids$MDD)), "coloc_pass MDD_DEG changed or is inconsistent")
assert_true(identical(as.integer(coloc_pass$BD_DEG), as.integer(coloc_pass$gene_id %in% disorder_ids$BD)), "coloc_pass BD_DEG changed or is inconsistent")
assert_true(nrow(coloc_pass) == 100L, "coloc_pass does not contain 100 rows")
assert_true(nrow(coloc_summary) == 24L, "coloc_summary does not contain 24 rows")

nominal_qc <- read_table(file.path(source_tables_dir, "nominal_BH05_qc.csv"))

excel_number_formats <- function(dt) {
  names_dt <- names(dt)
  formats <- setNames(rep(NA_character_, length(names_dt)), names_dt)
  lower <- tolower(names_dt)
  integer_cols <- grepl(
    paste(c(
      "^n_", "_n_", "^num_", "_count$", "count$", "rank$", "distance$",
      "samples$", "^nsnps$", "^cs95_n_snp$", "listed$", "rows$", "genes$",
      "variants$", "loci$", "tested$", "saved$", "input$", "error$", "null$",
      "workers$", "threads$", "pcs$", "members$"
    ), collapse = "|"),
    lower
  )
  flag_cols <- grepl(
    paste(c(
      "^deg$", "_deg$", "_eqtl$", "^lead_variant_in_", "^lead_variant_gwas_",
      "^gene_gwas_list$", "^gwas_support_"
    ), collapse = "|"),
    lower
  )
  scientific_cols <- lower %in% c("pp0", "pp1", "pp2") | grepl("pval|qval|fdr|gwasp$", lower)
  probability_cols <- lower %in% c(
    "pp3", "pp4", "pp34", "pp4_over_pp34", "lead_snp_pph4",
    "lead_snp_ppsho", "af", "case_fraction_s"
  ) | grepl("frac|fraction|si_min|impinfo", lower)
  decimal_cols <- grepl("slope|beta|shape|df$|se$|elapsed_sec|min_abs|dosage_r", lower)
  formats[integer_cols] <- "#,##0"
  formats[flag_cols] <- "0"
  formats[decimal_cols] <- "0.000"
  formats[probability_cols] <- "0.0000"
  formats[scientific_cols] <- "0.00E+00"
  formats[!vapply(dt, is.numeric, logical(1))] <- NA_character_
  formats
}

excel_widths <- function(dt, max_rows = 200L) {
  sample_count <- min(nrow(dt), max_rows)
  sample_dt <- if (sample_count > 0L) dt[seq_len(sample_count)] else dt
  widths <- vapply(names(dt), function(column) {
    values <- as.character(sample_dt[[column]])
    max(c(nchar(column), nchar(values)), na.rm = TRUE) + 2L
  }, numeric(1))
  widths <- pmin(pmax(widths, 8), 42)
  wide_cols <- grepl("genes$|gene_ids$|path|cache|thread_env|rule", names(dt), ignore.case = TRUE)
  widths[wide_cols] <- pmin(pmax(widths[wide_cols], 24), 80)
  widths
}

apply_excel_number_formats <- function(workbook, sheet, dt) {
  if (nrow(dt) == 0L || ncol(dt) == 0L) return(invisible(FALSE))
  formats <- excel_number_formats(dt)
  formats <- formats[!is.na(formats)]
  for (format in unique(unname(formats))) {
    columns <- match(names(formats)[formats == format], names(dt))
    style <- openxlsx::createStyle(numFmt = format)
    openxlsx::addStyle(
      workbook,
      sheet = sheet,
      style = style,
      rows = seq.int(2L, nrow(dt) + 1L),
      cols = columns,
      gridExpand = TRUE,
      stack = TRUE
    )
  }
  invisible(TRUE)
}

remove_dangling_openxlsx_relationships <- function(out_file) {
  if (!requireNamespace("xml2", quietly = TRUE)) fail("xml2 is required for XLSX validation")
  if (!requireNamespace("zip", quietly = TRUE)) fail("zip is required for XLSX validation")
  unpack_dir <- tempfile("xlsx-unpack-")
  repaired_file <- tempfile(fileext = ".xlsx")
  dir.create(unpack_dir)
  on.exit(unlink(c(unpack_dir, repaired_file), recursive = TRUE, force = TRUE), add = TRUE)
  unzip(out_file, exdir = unpack_dir)
  rel_dir <- file.path(unpack_dir, "xl", "worksheets", "_rels")
  rel_files <- if (dir.exists(rel_dir)) list.files(rel_dir, pattern = "[.]rels$", full.names = TRUE) else character()
  for (rel_file in rel_files) {
    document <- xml2::read_xml(rel_file)
    rel_nodes <- xml2::xml_find_all(document, ".//*[local-name()='Relationship']")
    source_dir <- dirname(dirname(rel_file))
    for (node in rel_nodes) {
      target <- xml2::xml_attr(node, "Target")
      rel_type <- xml2::xml_attr(node, "Type")
      target_path <- normalizePath(file.path(source_dir, target), mustWork = FALSE)
      is_drawing <- grepl("/(drawing|vmlDrawing)$", rel_type)
      if (is_drawing && !file.exists(target_path)) xml2::xml_remove(node)
    }
    xml2::write_xml(document, rel_file, options = "format")
  }
  zip::zip(
    repaired_file,
    files = list.files(unpack_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE),
    root = unpack_dir,
    include_directories = FALSE
  )
  if (!file.copy(repaired_file, out_file, overwrite = TRUE)) fail("Failed to repair workbook: ", out_file)
  invisible(out_file)
}

write_final_workbook <- function(sheets, out_file) {
  workbook <- openxlsx::createWorkbook()
  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fgFill = "#D9EAF7",
    border = "Bottom",
    borderStyle = "thin"
  )
  for (sheet in names(sheets)) {
    dt <- as.data.table(sheets[[sheet]])
    openxlsx::addWorksheet(workbook, sheet, gridLines = TRUE)
    openxlsx::writeDataTable(
      workbook,
      sheet = sheet,
      x = dt,
      tableStyle = "TableStyleLight9",
      headerStyle = header_style,
      withFilter = TRUE
    )
    openxlsx::freezePane(workbook, sheet, firstRow = TRUE)
    openxlsx::setColWidths(workbook, sheet, cols = seq_along(dt), widths = excel_widths(dt))
    apply_excel_number_formats(workbook, sheet, dt)
  }
  openxlsx::saveWorkbook(workbook, out_file, overwrite = TRUE)
  remove_dangling_openxlsx_relationships(out_file)
  invisible(out_file)
}

write_nominal_workbook <- function(nominal_dt, summary_dt, qc_dt, out_file) {
  workbook <- openxlsx::createWorkbook()
  for (sheet in c("nominal_BH05", "summary", "qc")) openxlsx::addWorksheet(workbook, sheet)
  openxlsx::writeDataTable(workbook, "nominal_BH05", nominal_dt)
  openxlsx::writeDataTable(workbook, "summary", summary_dt)
  openxlsx::writeDataTable(workbook, "qc", qc_dt)
  for (sheet in c("nominal_BH05", "summary", "qc")) openxlsx::freezePane(workbook, sheet, firstRow = TRUE)
  openxlsx::saveWorkbook(workbook, out_file, overwrite = TRUE)
  remove_dangling_openxlsx_relationships(out_file)
  invisible(out_file)
}

excel_column_number <- function(label) {
  letters <- utf8ToInt(toupper(label)) - utf8ToInt("A") + 1L
  sum(letters * 26L ^ rev(seq_along(letters) - 1L))
}

verify_workbook_structure <- function(out_file, sheets) {
  actual_sheets <- openxlsx::getSheetNames(out_file)
  if (!identical(actual_sheets, names(sheets))) fail("Unexpected sheets in ", out_file)
  results <- list()
  for (sheet_index in seq_along(sheets)) {
    sheet_name <- names(sheets)[[sheet_index]]
    expected <- as.data.table(sheets[[sheet_name]])
    header <- openxlsx::read.xlsx(
      out_file,
      sheet = sheet_name,
      rows = 1L,
      colNames = FALSE,
      skipEmptyRows = FALSE,
      skipEmptyCols = FALSE
    )
    observed_columns <- as.character(unlist(header[1L, ], use.names = FALSE))
    if (!identical(observed_columns, names(expected))) fail("Unexpected columns in ", out_file, " / ", sheet_name)
    rel_member <- file.path("xl", "worksheets", "_rels", paste0("sheet", sheet_index, ".xml.rels"))
    rel_connection <- unz(out_file, rel_member, open = "rb")
    rel_text <- rawToChar(readBin(rel_connection, what = "raw", n = 8192L))
    close(rel_connection)
    target_match <- regexec("Target=\"[.][.]/tables/([^\"]+[.]xml)\"", rel_text, perl = TRUE)
    target_parts <- regmatches(rel_text, target_match)[[1]]
    if (length(target_parts) != 2L) fail("Cannot locate table XML in ", out_file, " / ", sheet_name)
    table_member <- file.path("xl", "tables", target_parts[[2]])
    connection <- unz(out_file, table_member, open = "rb")
    prefix <- rawToChar(readBin(connection, what = "raw", n = 8192L))
    close(connection)
    dimension_match <- regexec("<table[^>]+ref=\"([^\"]+)\"", prefix, perl = TRUE)
    dimension_parts <- regmatches(prefix, dimension_match)[[1]]
    if (length(dimension_parts) != 2L) fail("Cannot read used range in ", out_file, " / ", sheet_name)
    dimension <- dimension_parts[[2]]
    last_cell <- tail(strsplit(dimension, ":", fixed = TRUE)[[1]], 1L)
    last_column <- sub("[0-9]+$", "", last_cell)
    last_row <- as.integer(sub("^[A-Z]+", "", last_cell))
    if (last_row != nrow(expected) + 1L || excel_column_number(last_column) != ncol(expected)) {
      fail("Unexpected used range in ", out_file, " / ", sheet_name, ": ", dimension)
    }
    results[[sheet_name]] <- data.table(
      workbook = basename(out_file),
      sheet = sheet_name,
      rows = nrow(expected),
      columns = ncol(expected),
      used_range = dimension
    )
  }
  rbindlist(results)
}

nominal_workbook_sheets <- list(
  nominal_BH05 = new_tables[["nominal_BH05.csv.gz"]],
  summary = new_summaries[["nominal_BH05_summary.csv"]],
  qc = nominal_qc
)
eqtl_workbook_sheets <- list(
  cis_independent = cis_public,
  nominal_BH05 = nominal_public,
  eQTL_summary = eqtl_summary_public,
  intersecting_genes = intersections_public
)
coloc_workbook_sheets <- list(
  coloc_pass = coloc_pass,
  coloc_summary = coloc_summary
)

candidate_nominal_xlsx <- file.path(candidate_tables_dir, "nominal_BH05.xlsx")
candidate_eqtl_xlsx <- file.path(candidate_final_dir, "eqtl_results.xlsx")
candidate_coloc_xlsx <- file.path(candidate_final_dir, "coloc_results.xlsx")
write_nominal_workbook(
  nominal_workbook_sheets$nominal_BH05,
  nominal_workbook_sheets$summary,
  nominal_workbook_sheets$qc,
  candidate_nominal_xlsx
)
write_final_workbook(eqtl_workbook_sheets, candidate_eqtl_xlsx)
write_final_workbook(coloc_workbook_sheets, candidate_coloc_xlsx)
workbook_validation <- rbindlist(list(
  verify_workbook_structure(candidate_nominal_xlsx, nominal_workbook_sheets),
  verify_workbook_structure(candidate_eqtl_xlsx, eqtl_workbook_sheets),
  verify_workbook_structure(candidate_coloc_xlsx, coloc_workbook_sheets)
))
fwrite(workbook_validation, file.path(report_dir, "workbook_validation.csv"))

impact_table <- rbindlist(impact_tables, use.names = TRUE, fill = TRUE)
impact_context <- rbindlist(impact_contexts, use.names = TRUE, fill = TRUE)
setcolorder(impact_context, c(
  "table", "context", "rows", "old_deg_rows", "new_deg_rows", "removed_rows",
  "old_deg_genes", "new_deg_genes", "removed_genes"
))
impact_overlap <- rbindlist(impact_overlaps, use.names = TRUE, fill = TRUE)
impact_overlap[, `:=`(
  removed_deg_gwas_rows = old_deg_gwas_rows - new_deg_gwas_rows,
  removed_deg_gwas_genes = old_deg_gwas_genes - new_deg_gwas_genes
)]
summary_impact <- rbindlist(summary_impacts, use.names = TRUE, fill = TRUE)
removed_gene_table <- as.data.table(legacy_degs$global)[gene_id %in% setdiff(legacy_ids, authoritative_ids)]

fwrite(impact_table, file.path(report_dir, "impact_by_table.csv"))
fwrite(impact_context, file.path(report_dir, "impact_by_context.csv"))
fwrite(impact_overlap, file.path(report_dir, "impact_by_disorder_overlap.csv"))
fwrite(summary_impact, file.path(report_dir, "impact_summary_counts.csv"))
fwrite(removed_gene_table, file.path(report_dir, "removed_legacy_DEGs.csv"))

validation <- data.table(
  metric = c(
    "sig_df_rows", "authoritative_gene_ids", "authoritative_gene_names", "legacy_gene_ids",
    "removed_legacy_genes", "added_genes", "MDD_DEG_genes", "BD_DEG_genes",
    "canonical_pair_rows", "canonical_signals", "canonical_index_rows", "canonical_recovered_ties",
    "boxplot_deg_pair_rows", "coloc_pass_rows", "coloc_summary_rows", "pre_change_snapshot_files",
    "workspace_files_matching_snapshot", "workspace_files_matching_candidate"
  ),
  value = c(
    nrow(sig_df), length(authoritative_ids), length(authoritative_names), length(legacy_ids),
    length(setdiff(legacy_ids, authoritative_ids)), length(setdiff(authoritative_ids, legacy_ids)),
    length(disorder_ids$MDD), length(disorder_ids$BD), nrow(canonical), uniqueN(canonical$signal_id),
    canonical[is_index_variant == TRUE, .N], canonical[recovery_status == "recovered_exact_tie", .N],
    nrow(boxplot_new), nrow(coloc_pass), nrow(coloc_summary), length(snapshot_targets),
    workspace_snapshot_matches, workspace_candidate_matches
  )
)
fwrite(validation, file.path(report_dir, "validation_report.csv"))

install_map <- c(
  setNames(file.path(tables_dir, c(row_files, summary_files, helper_files, slide_file)),
           file.path(candidate_tables_dir, c(row_files, summary_files, helper_files, slide_file))),
  setNames(file.path(tables_dir, "nominal_BH05.xlsx"), file.path(candidate_tables_dir, "nominal_BH05.xlsx")),
  setNames(file.path(final_dir, c("eqtl_results.xlsx", "coloc_results.xlsx")),
           file.path(candidate_final_dir, c("eqtl_results.xlsx", "coloc_results.xlsx")))
)

if (install_outputs) {
  for (candidate in names(install_map)) atomic_replace(candidate, install_map[[candidate]])
  cat("Installed validated annotation-only DEG outputs.\n")
} else {
  cat("Staged and validated annotation-only DEG outputs at:\n", output_dir, "\n")
  cat("Inspect the staged workbook previews before using --install.\n")
}
