#!/usr/bin/env Rscript

## compare frozen no-23andMe outputs with integrated MDD/BD GWAS outputs.
suppressPackageStartupMessages({
  library(data.table)
  library(here)
})

here::i_am(".git/HEAD")
repo_root <- here()
source(here("code", "11_eQTL_coloc", "utils.R"), chdir = FALSE)

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
  hit <- which(args == flag)
  if (!length(hit) || hit[[1]] == length(args)) return(default)
  args[[hit[[1]] + 1L]]
}

run_root <- here("processed-data", "11_eQTL_coloc", "seurat")
archive_root <- arg_value(
  "--archive",
  file.path(run_root, "archive", "no23andMe_2026-07-28")
)
report_file <- arg_value(
  "--report",
  here("code", "11_eQTL_coloc", "GWAS-23andMe-MDD-BD-impact-on-tables.md")
)
comparison_dir <- arg_value(
  "--comparison-dir",
  file.path(run_root, "comparison", "23andMe_2026-07-28")
)

if (!dir.exists(archive_root)) stop("Missing no-23andMe archive: ", archive_root)
dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)

## formatting helpers keep generated Markdown deterministic and ASCII-only.
fmt_int <- function(x) format(as.integer(x), big.mark = ",", scientific = FALSE, trim = TRUE)
fmt_num <- function(x, digits = 3L) {
  ifelse(is.na(x), "NA", format(round(x, digits), nsmall = digits, trim = TRUE))
}
fmt_delta <- function(x) ifelse(x > 0, paste0("+", fmt_int(x)), fmt_int(x))
collapse_values <- function(x) {
  x <- sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))])))
  if (length(x)) paste(x, collapse = ", ") else ""
}
md_escape <- function(x) gsub("|", "\\|", as.character(x), fixed = TRUE)

md_table <- function(dt) {
  dt <- as.data.table(dt)
  if (!nrow(dt)) return("(none)")
  vals <- lapply(dt, function(x) md_escape(ifelse(is.na(x), "NA", x)))
  rows <- do.call(paste, c(vals, sep = " | "))
  c(
    paste0("| ", paste(names(dt), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(dt)), collapse = " | "), " |"),
    paste0("| ", rows, " |")
  )
}

read_required <- function(path) {
  if (!file.exists(path)) stop("Missing required comparison input: ", path)
  fread(path)
}

write_comparison <- function(dt, name) {
  fwrite(dt, file.path(comparison_dir, name), sep = "\t", quote = FALSE, na = "NA")
  invisible(dt)
}

## compare stable eQTL columns after excluding GWAS annotations only.
core_comparison <- function(old, new, table_name, key_cols) {
  common <- intersect(names(old), names(new))
  core_cols <- common[!grepl("gwas|GWAS", common)]
  key_cols <- intersect(key_cols, core_cols)
  if (!length(key_cols)) stop("No comparison keys for ", table_name)

  old_core <- copy(old[, ..core_cols])
  new_core <- copy(new[, ..core_cols])
  setorderv(old_core, key_cols)
  setorderv(new_core, key_cols)
  comparison_tolerance <- 1e-12
  equal <- isTRUE(all.equal(
    old_core,
    new_core,
    check.attributes = FALSE,
    tolerance = comparison_tolerance
  ))
  data.table(
    table = table_name,
    old_rows = nrow(old),
    new_rows = nrow(new),
    row_delta = nrow(new) - nrow(old),
    stable_eqtl_deg_columns_equal = equal,
    numeric_tolerance = comparison_tolerance
  )
}

## count exact-variant and variant-or-curated-gene annotations by DEG scope.
overlap_metrics <- function(dt, version, table_name) {
  out <- list()
  base_scopes <- list(all = rep(TRUE, nrow(dt)), broad_DEG = dt$DEG == 1L)
  for (dis in c("MDD", "BD")) {
    scopes <- base_scopes
    dis_deg_col <- paste0(dis, "_DEG")
    if (dis_deg_col %in% names(dt)) scopes[[paste0(dis, "_DEG")]] <- dt[[dis_deg_col]] == 1L
    levels <- c(strict = "strict", suggestive_p1e5 = "exp")
    for (level_name in names(levels)) {
      suffix <- levels[[level_name]]
      flags <- c(
        exact_variant = paste0(dis, "_gwasVar_", suffix),
        variant_or_curated_gene = paste0(dis, "_gwas_", suffix)
      )
      for (annotation in names(flags)) {
        flag <- flags[[annotation]]
        if (!flag %in% names(dt)) next
        for (scope_name in names(scopes)) {
          keep <- !is.na(scopes[[scope_name]]) & scopes[[scope_name]] & dt[[flag]] == 1L
          out[[length(out) + 1L]] <- data.table(
            version = version,
            table = table_name,
            disorder = dis,
            level = level_name,
            annotation = annotation,
            scope = scope_name,
            n_rows = sum(keep),
            n_variants = uniqueN(dt$variant_id[keep]),
            n_genes = uniqueN(dt$gene_id[keep]),
            genes = collapse_values(dt$gene_name[keep])
          )
        }
      }
    }
  }
  rbindlist(out, use.names = TRUE, fill = TRUE)
}

## list exact variant annotations gained or lost on otherwise invariant rows.
annotation_changes <- function(old, new, table_name) {
  keys <- intersect(c("dataset_id", "context", "split", "gene_id", "gene_name", "variant_id"), names(old))
  out <- list()
  for (dis in c("MDD", "BD")) {
    levels <- c(strict = "strict", suggestive_p1e5 = "exp")
    for (level_name in names(levels)) {
      suffix <- levels[[level_name]]
      flag <- paste0(dis, "_gwasVar_", suffix)
      if (!flag %in% names(old) || !flag %in% names(new)) next
      a <- unique(old[get(flag) == 1L, ..keys])
      b <- unique(new[get(flag) == 1L, ..keys])
      a[, old_hit := TRUE]
      b[, new_hit := TRUE]
      z <- merge(a, b, by = keys, all = TRUE)
      z[, status := fcase(
        is.na(old_hit), "gained",
        is.na(new_hit), "lost",
        default = "retained"
      )]
      z[, `:=`(table = table_name, disorder = dis, level = level_name)]
      out[[length(out) + 1L]] <- z[status != "retained"]
    }
  }
  rbindlist(out, use.names = TRUE, fill = TRUE)
}

## load old and new 03-series tables.
table_specs <- list(
  map_significant_pairs = list(
    rel = file.path("tables", "map_significant_pairs.csv.gz"),
    keys = c("source", "dataset_id", "context", "split", "gene_id", "variant_id")
  ),
  nominal_BH05 = list(
    rel = file.path("tables", "nominal_BH05.csv.gz"),
    keys = c("dataset_id", "context", "split", "gene_id", "variant_id")
  ),
  map_GWASx = list(
    rel = file.path("tables", "map_sig_pairs_GWAS_relaxed.csv.gz"),
    keys = c("source", "dataset_id", "context", "split", "gene_id", "variant_id")
  )
)

core_checks <- list()
metric_parts <- list()
change_parts <- list()
for (nm in names(table_specs)) {
  spec <- table_specs[[nm]]
  old <- read_required(file.path(archive_root, spec$rel))
  new <- read_required(file.path(run_root, spec$rel))

  ## legacy 03c output used BPD labels; normalize only those column names.
  if (identical(nm, "map_GWASx")) {
    setnames(old, names(old), sub("BPD", "BD", names(old), fixed = TRUE))
    old <- rename_eqtl_overlap_cols(old)
    new <- rename_eqtl_overlap_cols(new)
  }

  core_checks[[nm]] <- core_comparison(old, new, nm, spec$keys)
  metric_parts[[paste0(nm, "_old")]] <- overlap_metrics(old, "no23andMe", nm)
  metric_parts[[paste0(nm, "_new")]] <- overlap_metrics(new, "full23andMe", nm)
  change_parts[[nm]] <- annotation_changes(old, new, nm)
}

core_checks <- rbindlist(core_checks)
overlap_long <- rbindlist(metric_parts, use.names = TRUE, fill = TRUE)
annotation_change_rows <- rbindlist(change_parts, use.names = TRUE, fill = TRUE)

overlap_wide <- dcast(
  overlap_long,
  table + disorder + level + annotation + scope ~ version,
  value.var = c("n_rows", "n_variants", "n_genes"),
  fill = 0L
)
overlap_genes_wide <- dcast(
  overlap_long,
  table + disorder + level + annotation + scope ~ version,
  value.var = "genes",
  fill = ""
)
setnames(
  overlap_genes_wide,
  c("no23andMe", "full23andMe"),
  c("genes_no23andMe", "genes_full23andMe")
)
overlap_wide <- merge(
  overlap_wide,
  overlap_genes_wide,
  by = c("table", "disorder", "level", "annotation", "scope"),
  all.x = TRUE,
  sort = FALSE
)
for (measure in c("n_rows", "n_variants", "n_genes")) {
  old_col <- paste0(measure, "_no23andMe")
  new_col <- paste0(measure, "_full23andMe")
  overlap_wide[, (paste0(measure, "_delta")) := get(new_col) - get(old_col)]
}

change_summary <- annotation_change_rows[, .(
  n_rows = .N,
  n_variants = uniqueN(variant_id),
  n_genes = uniqueN(gene_id),
  genes = collapse_values(gene_name)
), by = .(table, disorder, level, status)]

write_comparison(core_checks, "03_core_invariance.tsv")
write_comparison(overlap_long, "03_overlap_metrics_long.tsv")
write_comparison(overlap_wide, "03_overlap_metrics_comparison.tsv")
write_comparison(annotation_change_rows, "03_exact_variant_annotation_changes.tsv")
write_comparison(change_summary, "03_exact_variant_annotation_change_summary.tsv")

## verify that curated gene-list support did not change between runs.
file_sha256 <- function(path) {
  out <- system2("sha256sum", shQuote(path), stdout = TRUE)
  strsplit(out[[1]], "[[:space:]]+")[[1]][[1]]
}
gene_list_specs <- CJ(disorder = c("BD", "MDD"), list = c("gene_list", "prio_gene_list"))
gene_list_checks <- gene_list_specs[, {
  filename <- sprintf("GWAS_%s_%s.tsv", disorder, list)
  old_path <- file.path(archive_root, "gene_lists", filename)
  new_path <- here("processed-data", "ref", "GWAS", disorder, filename)
  old <- read_required(old_path)
  new <- read_required(new_path)
  .(
    n_genes_no23 = uniqueN(old$gene_symbol),
    n_genes_full = uniqueN(new$gene_symbol),
    sha256_identical = identical(file_sha256(old_path), file_sha256(new_path))
  )
}, by = .(disorder, list)]
write_comparison(gene_list_checks, "curated_GWAS_gene_list_invariance.tsv")

## compare sparse matched target-GWAS caches at a common p <= 1e-5 cutoff.
sparse_files <- data.table(
  disorder = c("BD", "MDD", "BD", "MDD"),
  version = rep(c("no23andMe", "full23andMe"), each = 2L),
  path = c(
    file.path(archive_root, "gwas_caches", "ref_GWAS", "GWAS-BD_flt_p1e-5_SI0.8_merged_maf05_tqtl-matched.tab.gz"),
    file.path(archive_root, "gwas_caches", "ref_GWAS", "GWAS-MDD_flt_p1e-5_SI0.8_merged_maf05_tqtl-matched.tab.gz"),
    here("processed-data", "ref", "GWAS", "GWAS-BD_full23andMe_preDENTIST_flt_p1e-5_SInone_merged_maf05_tqtl-matched.tab.gz"),
    here("processed-data", "ref", "GWAS", "GWAS-MDD_full23andMe_flt_p1e-5_SInone_merged_maf05_tqtl-matched.tab.gz")
  )
)
missing_sparse <- sparse_files[!file.exists(path)]
if (nrow(missing_sparse)) {
  stop("Missing common-threshold sparse cache(s): ", paste(missing_sparse$path, collapse = ", "))
}
sparse_counts <- sparse_files[, {
  x <- fread(path)
  .(
    matched_p1e5 = nrow(x),
    matched_p5e8 = sum(is.finite(x$p) & x$p <= GWAS_STRICT_P_THRESHOLD),
    exact_matches = sum(x$match_mode == "exact", na.rm = TRUE),
    swapped_matches = sum(x$match_mode == "swapped", na.rm = TRUE)
  )
}, by = .(disorder, version)]
sparse_wide <- dcast(sparse_counts, disorder ~ version, value.var = c("matched_p1e5", "matched_p5e8"))
sparse_wide[, `:=`(
  matched_p1e5_delta = matched_p1e5_full23andMe - matched_p1e5_no23andMe,
  matched_p5e8_delta = matched_p5e8_full23andMe - matched_p5e8_no23andMe
)]
write_comparison(sparse_counts, "sparse_matched_GWAS_counts.tsv")

## compare 05 coloc tables, including sensitivity-gated strong calls.
old_coloc <- read_required(file.path(archive_root, "coloc", "tables", "coloc_abf_results_gated.tsv.gz"))
new_coloc <- read_required(file.path(run_root, "coloc", "tables", "coloc_abf_results_gated.tsv.gz"))
old_coloc <- old_coloc[disorder %in% c("MDD", "BD")]
new_coloc <- new_coloc[disorder %in% c("MDD", "BD")]
coloc_key <- c("disorder", "dataset_id", "context", "gene_id", "gene_name")

coloc_summary_one <- function(dt, version) {
  dt[, .(
    version = version,
    n_candidate_loci = .N,
    n_candidate_genes = uniqueN(gene_id),
    n_raw_strong = sum(cat_raw == "strong_coloc", na.rm = TRUE),
    n_gated_strong = sum(cat_gated == "strong_coloc", na.rm = TRUE),
    gated_strong_genes = collapse_values(gene_name[cat_gated == "strong_coloc"])
  ), by = .(disorder, dataset_id, context)]
}
coloc_context <- rbindlist(list(
  coloc_summary_one(old_coloc, "no23andMe"),
  coloc_summary_one(new_coloc, "full23andMe")
))

coloc_context_wide <- dcast(
  coloc_context,
  disorder + dataset_id + context ~ version,
  value.var = c("n_candidate_loci", "n_raw_strong", "n_gated_strong"),
  fill = 0L
)
for (measure in c("n_candidate_loci", "n_raw_strong", "n_gated_strong")) {
  coloc_context_wide[, (paste0(measure, "_delta")) :=
    get(paste0(measure, "_full23andMe")) - get(paste0(measure, "_no23andMe"))]
}

coloc_total <- coloc_context[, .(
  n_candidate_loci = sum(n_candidate_loci),
  n_raw_strong = sum(n_raw_strong),
  n_gated_strong = sum(n_gated_strong),
  gated_strong_genes = collapse_values(unlist(strsplit(gated_strong_genes, ", ", fixed = TRUE)))
), by = .(disorder, version)]
coloc_total_wide <- dcast(
  coloc_total,
  disorder ~ version,
  value.var = c("n_candidate_loci", "n_raw_strong", "n_gated_strong")
)
for (measure in c("n_candidate_loci", "n_raw_strong", "n_gated_strong")) {
  coloc_total_wide[, (paste0(measure, "_delta")) :=
    get(paste0(measure, "_full23andMe")) - get(paste0(measure, "_no23andMe"))]
}

strong_sets <- function(dt, version) {
  out <- unique(dt[cat_gated == "strong_coloc", ..coloc_key])
  out[, (version) := TRUE]
  out
}
old_strong <- strong_sets(old_coloc, "no23andMe")
new_strong <- strong_sets(new_coloc, "full23andMe")
strong_changes <- merge(old_strong, new_strong, by = coloc_key, all = TRUE)
strong_changes[, status := fcase(
  is.na(no23andMe), "gained",
  is.na(full23andMe), "lost",
  default = "retained"
)]
strong_changes <- strong_changes[status != "retained"]
strong_change_summary <- strong_changes[, .(
  n_loci = .N,
  n_genes = uniqueN(gene_id),
  genes = collapse_values(gene_name)
), by = .(disorder, status)]

## final report sets contain only strong coloc results that passed sensitivity.
final_coloc_sets <- rbindlist(list(
  old_coloc[cat_gated == "strong_coloc", .(
    version = "no23",
    n_results = .N,
    n_genes = uniqueN(gene_id),
    genes = collapse_values(gene_name)
  ), by = disorder],
  new_coloc[cat_gated == "strong_coloc", .(
    version = "full",
    n_results = .N,
    n_genes = uniqueN(gene_id),
    genes = collapse_values(gene_name)
  ), by = disorder]
))
final_coloc_wide <- dcast(
  final_coloc_sets,
  disorder ~ version,
  value.var = c("n_results", "n_genes"),
  fill = 0L
)
final_coloc_wide[, `:=`(
  result_delta = n_results_full - n_results_no23,
  gene_delta = n_genes_full - n_genes_no23
)]

## quantify posterior changes only for loci present in both candidate tables.
shared <- merge(
  old_coloc[, c(coloc_key, "PP4", "cat_raw", "cat_gated", "lead_snp"), with = FALSE],
  new_coloc[, c(coloc_key, "PP4", "cat_raw", "cat_gated", "lead_snp"), with = FALSE],
  by = coloc_key,
  suffixes = c("_no23", "_full")
)
shared_pp <- shared[, .(
  n_shared_candidates = .N,
  pearson_PP4 = if (.N > 1L) cor(PP4_no23, PP4_full, use = "complete.obs") else NA_real_,
  median_abs_PP4_change = median(abs(PP4_full - PP4_no23), na.rm = TRUE),
  max_abs_PP4_change = max(abs(PP4_full - PP4_no23), na.rm = TRUE),
  lead_snp_unchanged_pct = 100 * mean(lead_snp_no23 == lead_snp_full, na.rm = TRUE)
), by = disorder]

category_transitions <- shared[, .N, by = .(
  disorder,
  cat_no23 = fifelse(is.na(cat_raw_no23), "none", cat_raw_no23),
  cat_full = fifelse(is.na(cat_raw_full), "none", cat_raw_full)
)][order(disorder, -N, cat_no23, cat_full)]

write_comparison(coloc_context, "coloc_context_long.tsv")
write_comparison(coloc_context_wide, "coloc_context_comparison.tsv")
write_comparison(coloc_total, "coloc_disorder_totals_long.tsv")
write_comparison(strong_changes, "coloc_gated_strong_changes.tsv")
write_comparison(strong_change_summary, "coloc_gated_strong_change_summary.tsv")
write_comparison(shared_pp, "coloc_shared_candidate_PP4_changes.tsv")
write_comparison(category_transitions, "coloc_category_transitions.tsv")

## collect release provenance from current and archived run metadata.
meta_files <- rbindlist(lapply(c("MDD", "BD"), function(dis) {
  data.table(
    disorder = dis,
    version = c("no23andMe", "full23andMe"),
    path = c(
      file.path(archive_root, "coloc", dis, "coloc_astro.runmeta.tsv.gz"),
      file.path(run_root, "coloc", dis, "coloc_astro.runmeta.tsv.gz")
    )
  )
}))
archive_sha <- fread(
  file.path(archive_root, "gwas_input_sha256.tsv"),
  header = FALSE,
  col.names = c("sha256", "path")
)
provenance <- meta_files[, {
  x <- fread(path)
  if ("gwas_bcf" %in% names(x)) {
    bcf_path <- x$gwas_bcf[[1]]
    bcf_sha <- x$gwas_bcf_sha256[[1]]
  } else {
    pattern <- if (disorder == "BD") "bip2024_eur_no23andMe.hg38.bcf" else "pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.bcf"
    old_input <- archive_sha[endsWith(path, pattern)]
    if (nrow(old_input) != 1L) stop("Cannot resolve archived GWAS provenance for ", disorder)
    bcf_path <- old_input$path[[1]]
    bcf_sha <- old_input$sha256[[1]]
  }
  .(
    release = if ("gwas_release" %in% names(x)) x$gwas_release[[1]] else "public_no23andMe",
    bcf = basename(bcf_path),
    sha256 = bcf_sha,
    SI_filter = if (is.na(x$gwas_si_min[[1]])) "none" else paste0(">=", x$gwas_si_min[[1]])
  )
}, by = .(disorder, version)]
write_comparison(provenance, "GWAS_input_provenance.tsv")

## build concise report-facing tables.
core_report <- copy(core_checks)
core_report[, `:=`(
  old_rows = fmt_int(old_rows),
  new_rows = fmt_int(new_rows),
  row_delta = vapply(row_delta, fmt_delta, character(1)),
  stable_eqtl_deg_columns_equal = ifelse(stable_eqtl_deg_columns_equal, "yes", "NO")
)]

sparse_report <- sparse_wide[, .(
  disorder,
  p1e5_no23 = fmt_int(matched_p1e5_no23andMe),
  p1e5_full = fmt_int(matched_p1e5_full23andMe),
  p1e5_delta = vapply(matched_p1e5_delta, fmt_delta, character(1)),
  p5e8_no23 = fmt_int(matched_p5e8_no23andMe),
  p5e8_full = fmt_int(matched_p5e8_full23andMe),
  p5e8_delta = vapply(matched_p5e8_delta, fmt_delta, character(1))
)]

overlap_report <- overlap_wide[
  table %in% c("map_significant_pairs", "nominal_BH05") &
    annotation == "exact_variant" & scope == "all" &
    (table == "map_significant_pairs" | level == "strict"),
  .(
    table, disorder, level,
    rows_no23 = fmt_int(n_rows_no23andMe),
    rows_full = fmt_int(n_rows_full23andMe),
    row_delta = vapply(n_rows_delta, fmt_delta, character(1)),
    genes_no23 = fmt_int(n_genes_no23andMe),
    genes_full = fmt_int(n_genes_full23andMe),
    gene_delta = vapply(n_genes_delta, fmt_delta, character(1))
  )
]

deg_overlap_report <- overlap_wide[
  table == "map_significant_pairs" &
    ((annotation == "exact_variant" & scope == "broad_DEG") |
      (annotation == "variant_or_curated_gene" & scope == paste0(disorder, "_DEG"))),
  .(
    disorder, level, annotation, scope,
    rows_no23 = fmt_int(n_rows_no23andMe),
    rows_full = fmt_int(n_rows_full23andMe),
    row_delta = vapply(n_rows_delta, fmt_delta, character(1)),
    gene_names_no23 = genes_no23andMe,
    gene_names_full = genes_full23andMe
  )
]

coloc_report <- final_coloc_wide[, .(
  disorder,
  final_results_no23 = fmt_int(n_results_no23),
  final_results_full = fmt_int(n_results_full),
  result_delta = vapply(result_delta, fmt_delta, character(1)),
  distinct_genes_no23 = fmt_int(n_genes_no23),
  distinct_genes_full = fmt_int(n_genes_full),
  gene_delta = vapply(gene_delta, fmt_delta, character(1))
)]

short_final_gene_sets_report <- final_coloc_sets[n_genes < 20L, .(
  disorder,
  version,
  final_results = fmt_int(n_results),
  distinct_genes = fmt_int(n_genes),
  gene_symbols = genes
)]

coloc_context_report <- coloc_context_wide[, .(
  disorder, context,
  final_no23 = fmt_int(n_gated_strong_no23andMe),
  final_full = fmt_int(n_gated_strong_full23andMe),
  delta = vapply(n_gated_strong_delta, fmt_delta, character(1))
)][delta != "0"]

strong_change_report <- strong_change_summary[, .(
  disorder, status,
  results = fmt_int(n_loci),
  distinct_genes = fmt_int(n_genes),
  gene_symbols = genes
)]

## report caveats are factual properties of the supplied integrated releases.
all_core_identical <- all(core_checks$stable_eqtl_deg_columns_equal)
report <- c(
  "# Impact of 23andMe-inclusive MDD and BD GWAS on eQTL/DEG tables",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Scope",
  "",
  "This comparison holds the 119-donor genotype data, tensorQTL eQTL results, DEG definitions, curated GWAS gene lists, coloc priors, and sensitivity criteria fixed. It changes only the MDD and BD GWAS inputs from public European no-23andMe statistics to reconstructed European statistics that include 23andMe.",
  "",
  "The deltas therefore measure the effect of switching supplied GWAS files, not an isolated marginal effect of the 23andMe cohorts. The integrated meta-analysis statistics and available row-level QC fields also differ from the public files; in particular, integrated MDD/BD has no combined imputation-quality (`SI`) field.",
  "",
  "## Terminology",
  "",
  "- `no23` means the archived public European GWAS that excludes 23andMe. `full` means the local 23andMe-inclusive European reconstruction obtained by meta-analyzing that public result with the delivered 23andMe component; it does not mean an official final paper file. MDD is a release-matched reconstruction. BD uses v7.0 associations with v7.2 annotations and is `pre-DENTIST`, meaning the paper's final post-meta-analysis LD-based DENTIST QC was not reproduced.",
  "- A `canonical variant` is identified as GRCh38 `chromosome:position:REF:ALT`. An `eQTL row` or `eQTL pair` is one tensorQTL association between such a variant and an eGene, meaning a gene whose expression is associated with that variant, in a cell context and analysis split. The eQTL result, not the GWAS, supplies the gene assignment.",
  "- `strict` or `exact_variant` means that the same canonical eQTL variant has GWAS P <= 5e-8. `suggestive_p1e5` means exact-variant GWAS P < 1e-5.",
  "- `broad_DEG` is membership in the project-wide union of the author-recommended DEG tables. `MDD_DEG` and `BD_DEG` are the diagnosis-specific NTC-versus-MDD and NTC-versus-BD DEG subsets.",
  "- A `curated GWAS gene` comes from a paper-derived GWAS gene list held fixed between runs. `variant_or_curated_gene` means exact variant support OR curated-gene-list membership, so it is not necessarily a variant-level overlap.",
  "",
  "## GWAS inputs",
  "",
  md_table(provenance),
  "",
  "The integrated MDD file is the 23andMe-inclusive European meta-analysis. The integrated BD file is the reconstructed 23andMe-inclusive European meta-analysis and remains pre-DENTIST; it is not represented here as the exact final paper release. SI is absent from both integrated BCFs, so no post-integration SI filter was applied. The archived public files used SI >= 0.8.",
  "",
  "## Cache and table validation",
  "",
  md_table(core_report),
  "",
  paste0("Stable non-GWAS eQTL and DEG columns were equal across all compared 03-series row tables at the documented 1e-12 numeric tolerance: ", ifelse(all_core_identical, "yes", "NO"), "."),
  "",
  "Curated broad and prioritized GWAS gene-list files were byte-identical between runs:",
  "",
  md_table(gene_list_checks),
  "",
  "Matched variants in the fixed 119-donor PLINK2 target:",
  "",
  md_table(sparse_report),
  "",
  "These are matched variant rows, not independent GWAS loci.",
  "",
  "## eQTL/GWAS and DEG-scoped overlap changes",
  "",
  "### Exact eQTL-GWAS variant overlaps (no DEG filter)",
  "",
  "This table counts eQTL rows whose exact canonical variant passes the stated GWAS threshold; `rows_*` are context-specific eQTL association rows and `genes_*` are distinct tensorQTL eGenes, with no DEG or curated-GWAS-gene requirement. `map_significant_pairs` is the significant cis/independent eQTL set, whereas `nominal_BH05` is the broader nominal eQTL set passing BH FDR 0.05.",
  "",
  md_table(overlap_report),
  "",
  "### DEG-scoped significant-eQTL overlaps",
  "",
  "This table restricts significant eQTL rows by DEG status: `exact_variant` with `broad_DEG` is the exact GWAS-variant/eQTL/DEG trifecta, while `variant_or_curated_gene` with a diagnosis-specific DEG permits either exact-variant support or curated-GWAS-gene support. `gene_names_*` lists the distinct eGenes in each intersection.",
  "",
  md_table(deg_overlap_report),
  "",
  "Exact gained/lost row and gene lists are in `processed-data/11_eQTL_coloc/seurat/comparison/23andMe_2026-07-28/03_exact_variant_annotation_changes.tsv`. Combined variant-or-curated-gene and DEG-scope counts are in `03_overlap_metrics_comparison.tsv` in the same directory.",
  "",
  "## Final colocalization changes",
  "",
  "Only final strong-coloc results that passed sensitivity testing are shown. One result is one disorder-by-cell-context-by-eGene combination, so result counts can exceed distinct-gene counts.",
  "",
  md_table(coloc_report),
  "",
  "Complete gene-symbol lists are shown below for final sets containing fewer than 20 distinct genes:",
  "",
  md_table(short_final_gene_sets_report),
  "",
  "Cell contexts with a changed final strong-coloc result count:",
  "",
  md_table(coloc_context_report),
  "",
  "Final strong-coloc gains and losses are shown with complete gene-symbol lists because every changed set contains fewer than 20 distinct genes:",
  "",
  md_table(strong_change_report),
  "",
  "## Reproducibility",
  "",
  paste0("Frozen archive: `", normalizePath(archive_root), "`"),
  "",
  paste0("Machine-readable comparison tables: `", normalizePath(comparison_dir), "`"),
  "",
  "Regenerate this report after rerunning 03-05 with:",
  "",
  "```bash",
  "Rscript code/11_eQTL_coloc/07_compare_23andme_impact.R",
  "```"
)

dir.create(dirname(report_file), recursive = TRUE, showWarnings = FALSE)
writeLines(report, report_file, useBytes = TRUE)
message("Wrote comparison report: ", report_file)
