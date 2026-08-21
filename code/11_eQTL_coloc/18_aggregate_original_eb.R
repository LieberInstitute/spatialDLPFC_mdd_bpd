#!/usr/bin/env Rscript

## aggregate EB sensitivity for the original 100 coloc_pass targets.

pilot_lib <- paste0(
  "/home/gpertea/work/R/spatialDLPFC_mdd_bpd/processed-data/11_eQTL_coloc/",
  "seurat/coloc/susie_exploratory/Rlib"
)
.libPaths(c(pilot_lib, .libPaths()))
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: 18_aggregate_original_eb.R BASE")
base <- normalizePath(args[[1L]], mustWork = TRUE)
out_dir <- file.path(base, "aggregate_eb_20260820")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
manifest <- fread(file.path(base, "target_manifest.tsv"))
mapk_dir <- paste0(
  "/home/gpertea/work/R/spatialDLPFC_mdd_bpd/processed-data/11_eQTL_coloc/seurat/",
  "coloc/susie_next_no23andMe_20260820/mapk3_eb"
)
mapk_old <- paste0(
  "/home/gpertea/work/R/spatialDLPFC_mdd_bpd/processed-data/11_eQTL_coloc/seurat/",
  "coloc/susie_targeted_no23andMe/MAPK3_pilot_20260819"
)

canonical_cs <- function(fit) {
  cs <- fit$sets$cs
  if (is.null(cs) || !length(cs)) return(character())
  sort(vapply(cs, function(i) paste(sort(names(fit$pip)[i]), collapse = ","), character(1L)))
}

collect_cs <- function(fit, target_id, trait) {
  cs <- fit$sets$cs
  if (is.null(cs) || !length(cs)) return(data.table())
  rbindlist(lapply(seq_along(cs), function(i) {
    idx <- cs[[i]]
    data.table(
      target_id = target_id, trait = trait, cs = names(cs)[i],
      variant_id = names(fit$pip)[idx], pip = fit$pip[idx],
      cs_size = length(idx), purity_min_abs_corr = fit$sets$purity$min.abs.corr[i]
    )
  }))
}

collect_h4 <- function(coloc_result, target_id) {
  result <- as.data.table(coloc_result$results)
  summary <- as.data.table(coloc_result$summary)
  pp_cols <- grep("^SNP.PP.H4", names(result), value = TRUE)
  if (!nrow(summary) || !length(pp_cols)) return(data.table())
  rbindlist(lapply(seq_along(pp_cols), function(i) {
    x <- data.table(variant_id = result$snp, h4_conditional_variant_pp = result[[pp_cols[[i]]]])
    setorder(x, -h4_conditional_variant_pp)
    x[, `:=`(rank = .I, cumulative_pp = cumsum(h4_conditional_variant_pp))]
    keep_n <- max(10L, which(x$cumulative_pp >= 0.95)[1L], na.rm = TRUE)
    x <- x[seq_len(min(nrow(x), keep_n))]
    x[, `:=`(
      target_id = target_id, signal_pair = i,
      hit1 = summary$hit1[[i]], hit2 = summary$hit2[[i]],
      pair_H3 = summary$PP.H3.abf[[i]], pair_H4 = summary$PP.H4.abf[[i]]
    )]
    x
  }), fill = TRUE)
}

summarize_fit <- function(fit, target_id, trait) {
  d <- fit$R_finite_diagnostics
  get1 <- function(name) if (is.null(d[[name]]) || !length(d[[name]])) NA else d[[name]][1L]
  data.table(
    target_id = target_id, trait = trait, converged = isTRUE(fit$converged),
    niter = fit$niter, n_credible_sets = length(fit$sets$cs),
    max_pip = max(fit$pip), top_variant = names(fit$pip)[which.max(fit$pip)],
    lambda_bias = get1("lambda_bias"), B_corrected = get1("B_corrected"),
    Q_art = get1("Q_art"), sensitivity_flag = get1("R_sensitivity_flag"),
    reliability_flag = get1("R_reliability_flag")
  )
}

rows <- list(fits = list(), comparisons = list(), coloc = list(), cs = list(),
             h4_variants = list(), exclusions = list(), ld = list(), timings = list(), warnings = list(),
             metadata = list())

for (i in seq_len(nrow(manifest))) {
  target <- manifest[i]
  id <- target$target_id
  if (isTRUE(target$mapk3_pilot_complete)) {
    key <- paste(target$disorder, target$dataset_id, sep = "__")
    x <- readRDS(file.path(mapk_dir, "fits", paste0(key, ".rds")))
    fits <- rbindlist(list(
      summarize_fit(x$fit_eqtl_wrapper, id, "eqtl"),
      summarize_fit(x$fit_gwas_eb_wrapper, id, "gwas")
    ), fill = TRUE)
    comparisons <- fread(file.path(mapk_dir, "direct_wrapper_comparison.tsv"))[
      comparison_key == key
    ][, target_id := id]
    coloc <- as.data.table(x$coloc_eb$summary)
    if (!nrow(coloc)) coloc <- data.table(status = "no_credible_set_pair")
    audit <- fread(file.path(mapk_old, "diagnostics", paste0(key, "_variant_audit.tsv.gz")))
    exclusions <- audit[, .N, by = .(
      exclusion_reason = fifelse(is.na(exclusion_reason), "included", exclusion_reason)
    )]
    ld <- data.table()
    timings <- fread(file.path(mapk_dir, "timings.tsv"))[comparison_key == key]
    warnings <- fread(file.path(mapk_dir, "warnings.tsv"), sep = "\t")[
      grepl(key, warning, fixed = TRUE)
    ]
    metadata <- data.table(
      target_id = id, source = "MAPK3_EB_reused_verified_eqtl",
      n_variants = length(x$ids), peak_rss_kb = NA_real_
    )
    eqtl_fit <- x$fit_eqtl_wrapper
    gwas_fit <- x$fit_gwas_eb_wrapper
  } else {
    x <- readRDS(file.path(base, "results_eb", paste0(id, ".rds")))
    fits <- copy(x$fit_summary)
    rename <- c(
      mismatch_lambda_bias = "lambda_bias", mismatch_B_corrected = "B_corrected",
      mismatch_Q_art = "Q_art", finite_sensitivity_flag = "sensitivity_flag",
      finite_reliability_flag = "reliability_flag"
    )
    present <- names(rename)[names(rename) %in% names(fits)]
    setnames(fits, present, unname(rename[present]))
    baseline_fit <- readRDS(file.path(base, "results", paste0(id, ".rds")))
    comparisons <- rbindlist(list(
      copy(x$direct_wrapper_comparison),
      data.table(
        trait = "eqtl_reuse_check",
        max_abs_pip_diff = max(abs(x$fit_eqtl$pip - baseline_fit$fit_eqtl$pip)),
        max_abs_alpha_diff = max(abs(x$fit_eqtl$alpha - baseline_fit$fit_eqtl$alpha)),
        final_elbo_diff = tail(x$fit_eqtl$elbo, 1L) - tail(baseline_fit$fit_eqtl$elbo, 1L),
        cs_identical = identical(canonical_cs(x$fit_eqtl), canonical_cs(baseline_fit$fit_eqtl))
      )
    ), fill = TRUE)
    coloc <- as.data.table(x$coloc_primary$summary)
    if (!nrow(coloc)) coloc <- data.table(status = "no_credible_set_pair")
    exclusions <- copy(x$exclusions)
    ld <- copy(x$ld_diagnostics)
    timings <- copy(x$timings)
    warnings <- data.table(warning = x$warnings)
    metadata <- copy(x$metadata)
    eqtl_fit <- x$fit_eqtl
    gwas_fit <- x$fit_gwas
  }
  for (value in c("fits", "comparisons", "coloc", "exclusions", "ld", "timings", "warnings")) {
    object <- get(value)
    object[, target_id := id]
    rows[[value]][[id]] <- object
  }
  rows$metadata[[id]] <- metadata
  rows$h4_variants[[id]] <- collect_h4(if (isTRUE(target$mapk3_pilot_complete)) x$coloc_eb else x$coloc_primary, id)
  rows$cs[[paste0(id, "__eqtl")]] <- collect_cs(eqtl_fit, id, "eqtl")
  rows$cs[[paste0(id, "__gwas")]] <- collect_cs(gwas_fit, id, "gwas")
}

tables <- lapply(rows, rbindlist, fill = TRUE)

## collect GNU time records across the adaptive worker pools.
extract_resource <- function(id) {
  candidates <- Sys.glob(file.path(base, "logs_eb*", paste0(id, ".err")))
  lines <- unlist(lapply(candidates, readLines, warn = FALSE), use.names = FALSE)
  extract_num <- function(pattern) {
    value <- grep(pattern, lines, value = TRUE)
    if (!length(value)) return(NA_real_)
    text <- trimws(sub("^[^:]+:", "", tail(value, 1L)))
    as.numeric(sub("[^0-9.eE+-].*$", "", text))
  }
  data.table(
    target_id = id, user_seconds = extract_num("User time \\(seconds\\)"),
    system_seconds = extract_num("System time \\(seconds\\)"),
    cpu_percent = extract_num("Percent of CPU"),
    peak_rss_kb_external = extract_num("Maximum resident set size")
  )
}
resources <- rbindlist(lapply(
  manifest[mapk3_pilot_complete == FALSE, target_id], extract_resource
))
mapk_peak <- as.numeric(sub("^VmHWM:[[:space:]]*([0-9]+).*$", "\\1", grep(
  "^VmHWM:", readLines(file.path(mapk_dir, "process_memory.txt")), value = TRUE
)))
resources <- rbindlist(list(
  resources,
  manifest[mapk3_pilot_complete == TRUE, .(
    target_id, user_seconds = NA_real_, system_seconds = NA_real_,
    cpu_percent = NA_real_, peak_rss_kb_external = mapk_peak
  )]
), fill = TRUE)
tables$resources <- resources
for (name in names(tables)) {
  fwrite(tables[[name]], file.path(out_dir, paste0(name, ".tsv.gz")), sep = "\t")
}

## compare target-level conclusions with the immutable baseline aggregate.
coloc <- tables$coloc
pair_summary <- coloc[, .(
  eb_n_pairs = sum(!is.na(PP.H4.abf)),
  eb_max_H3 = if (any(!is.na(PP.H3.abf))) max(PP.H3.abf, na.rm = TRUE) else NA_real_,
  eb_max_H4 = if (any(!is.na(PP.H4.abf))) max(PP.H4.abf, na.rm = TRUE) else NA_real_,
  eb_n_H3_ge_0_8 = sum(PP.H3.abf >= 0.8, na.rm = TRUE),
  eb_n_H4_ge_0_8 = sum(PP.H4.abf >= 0.8, na.rm = TRUE)
), by = target_id]
fit_wide <- dcast(
  tables$fits[, .(target_id, trait, n_credible_sets, max_pip, top_variant)],
  target_id ~ trait, value.var = c("n_credible_sets", "max_pip", "top_variant")
)
baseline <- fread(file.path(base, "aggregate", "target_summary.tsv.gz"))
comparison <- merge(baseline, pair_summary, by = "target_id", all.x = TRUE)
comparison <- merge(comparison, fit_wide, by = "target_id", all.x = TRUE)
fwrite(comparison, file.path(out_dir, "baseline_vs_eb_target_summary.tsv.gz"), sep = "\t")

audit <- data.table(
  check = c("manifest_rows", "eb_target_rows", "trait_fits", "nonconverged",
            "direct_wrapper_nonzero", "sensitivity_flags", "reliability_flags"),
  value = c(
    nrow(manifest), uniqueN(tables$fits$target_id), nrow(tables$fits),
    sum(!tables$fits$converged),
    sum(tables$comparisons$max_abs_pip_diff != 0, na.rm = TRUE),
    sum(tables$fits$sensitivity_flag %in% TRUE, na.rm = TRUE),
    sum(tables$fits$reliability_flag %in% TRUE, na.rm = TRUE)
  )
)
fwrite(audit, file.path(out_dir, "completion_audit.tsv"), sep = "\t")
message("aggregated original 100-target EB sensitivity")
