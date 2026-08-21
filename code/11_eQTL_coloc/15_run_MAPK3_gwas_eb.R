#!/usr/bin/env Rscript

## refit only the MAPK3 GWAS models with empirical-Bayes LD mismatch correction.
## reuse the previously verified eQTL fits and matched coloc datasets unchanged.

pilot_lib <- paste0(
  "/home/gpertea/work/R/spatialDLPFC_mdd_bpd/processed-data/11_eQTL_coloc/",
  "seurat/coloc/susie_exploratory/Rlib"
)
.libPaths(c(pilot_lib, .libPaths()))

suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(susieR)
})

data.table::setDTthreads(1L)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) stop("usage: 15_run_MAPK3_gwas_eb.R [PILOT_DIR] [OUT_DIR]")
pilot_dir <- if (length(args) >= 1L) args[[1L]] else paste0(
  "/home/gpertea/work/R/spatialDLPFC_mdd_bpd/processed-data/11_eQTL_coloc/seurat/",
  "coloc/susie_targeted_no23andMe/MAPK3_pilot_20260819"
)
out_dir <- if (length(args) == 2L) args[[2L]] else paste0(
  "/home/gpertea/work/R/spatialDLPFC_mdd_bpd/processed-data/11_eQTL_coloc/seurat/",
  "coloc/susie_next_no23andMe_20260820/mapk3_eb"
)
pilot_dir <- normalizePath(pilot_dir, mustWork = TRUE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, mustWork = TRUE)
dir.create(file.path(out_dir, "fits"), showWarnings = FALSE)

## use exactly the pilot settings, adding only the documented EB mismatch mode.
fit_files <- sort(list.files(file.path(pilot_dir, "fits"), "\\.rds$", full.names = TRUE))
if (length(fit_files) != 24L) stop("Expected 24 completed MAPK3 pilot fits")
warnings_seen <- character()
timing_rows <- list()

timed <- function(key, step, expr) {
  gc()
  started <- Sys.time()
  tm <- system.time(value <- withCallingHandlers(
    force(expr),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, paste(key, step, conditionMessage(w), sep = ": "))
      invokeRestart("muffleWarning")
    }
  ))
  timing_rows[[length(timing_rows) + 1L]] <<- data.table(
    comparison_key = key, step = step,
    started = format(started, "%Y-%m-%d %H:%M:%S %z"),
    elapsed_seconds = unname(tm[["elapsed"]]),
    user_seconds = unname(tm[["user.self"]]),
    system_seconds = unname(tm[["sys.self"]])
  )
  value
}

canonical_cs <- function(fit) {
  cs <- fit$sets$cs
  if (is.null(cs) || !length(cs)) return(character())
  sort(vapply(cs, function(i) paste(sort(names(fit$pip)[i]), collapse = ","), character(1L)))
}

compare_fits <- function(a, b) {
  data.table(
    max_abs_pip_diff = max(abs(a$pip - b$pip)),
    max_abs_alpha_diff = max(abs(a$alpha - b$alpha)),
    final_elbo_diff = tail(a$elbo, 1L) - tail(b$elbo, 1L),
    niter_a = a$niter, niter_b = b$niter,
    converged_a = isTRUE(a$converged), converged_b = isTRUE(b$converged),
    cs_identical = identical(canonical_cs(a), canonical_cs(b)),
    n_cs_a = length(a$sets$cs), n_cs_b = length(b$sets$cs)
  )
}

summarize_fit <- function(key, mode, fit) {
  d <- fit$R_finite_diagnostics
  data.table(
    comparison_key = key, mode = mode, converged = isTRUE(fit$converged),
    niter = fit$niter, n_credible_sets = length(fit$sets$cs),
    max_pip = max(fit$pip), top_variant = names(fit$pip)[which.max(fit$pip)],
    final_elbo = tail(fit$elbo, 1L), lambda_bias = d$lambda_bias %||% NA_real_,
    B_corrected = d$B_corrected %||% NA_real_, Q_art = d$Q_art %||% NA_real_,
    finite_effective_rank = d$effective_rank %||% NA_real_,
    finite_rank_over_B = d$r_over_B %||% NA_real_,
    finite_max_penalty = if (is.null(d$per_variable_penalty)) NA_real_ else max(d$per_variable_penalty),
    finite_sensitivity_flag = d$R_sensitivity_flag %||% NA,
    finite_reliability_flag = d$R_reliability_flag %||% NA
  )
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
fit_rows <- list()
comparison_rows <- list()
coloc_rows <- list()

for (fit_file in fit_files) {
  key <- sub("\\.rds$", "", basename(fit_file))
  final_file <- file.path(out_dir, "fits", basename(fit_file))
  message(format(Sys.time()), " | starting ", key)
  old <- readRDS(fit_file)
  if (!identical(old$ids, old$D_gwas$snp)) stop(key, ": pilot GWAS variant order changed")
  if (!identical(names(old$fit_eqtl_wrapper$pip), old$ids)) stop(key, ": eQTL fit order changed")

  eb_args <- modifyList(old$susie_args_gwas, list(R_mismatch = "eb"))
  fit_direct <- timed(
    key, "susie_rss_gwas_eb",
    do.call(susieR::susie_rss, c(list(
      z = old$D_gwas$z, R = old$D_gwas$LD, n = old$D_gwas$N
    ), eb_args))
  )
  fit_wrapper <- timed(
    key, "runsusie_gwas_eb",
    do.call(coloc::runsusie, c(list(d = old$D_gwas), eb_args))
  )
  coloc_eb <- timed(
    key, "coloc_susie_eb",
    coloc::coloc.susie(old$fit_eqtl_wrapper, fit_wrapper,
                       p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  )

  fit_rows[[paste0(key, "__none")]] <- summarize_fit(key, "none", old$fit_gwas_wrapper)
  fit_rows[[paste0(key, "__eb")]] <- summarize_fit(key, "eb", fit_wrapper)
  comparison_rows[[key]] <- cbind(
    data.table(comparison_key = key, comparison = "gwas_eb_susie_rss_vs_runsusie"),
    compare_fits(fit_direct, fit_wrapper)
  )
  cs <- as.data.table(coloc_eb$summary)
  if (!nrow(cs)) cs <- data.table(status = "no_credible_set_pair")
  cs[, comparison_key := key]
  coloc_rows[[key]] <- cs

  saveRDS(list(
    source_fit = normalizePath(fit_file), ids = old$ids,
    fit_eqtl_wrapper = old$fit_eqtl_wrapper,
    fit_gwas_none = old$fit_gwas_wrapper,
    fit_gwas_eb_direct = fit_direct, fit_gwas_eb_wrapper = fit_wrapper,
    coloc_none = old$coloc_prefit_primary, coloc_eb = coloc_eb,
    susie_args_gwas_eb = eb_args
  ), final_file, compress = "gzip")
}

## write small audit tables; the fit objects remain in the isolated data directory.
fwrite(rbindlist(fit_rows, fill = TRUE), file.path(out_dir, "fit_summary.tsv"), sep = "\t")
fwrite(rbindlist(comparison_rows, fill = TRUE), file.path(out_dir, "direct_wrapper_comparison.tsv"), sep = "\t")
fwrite(rbindlist(coloc_rows, fill = TRUE), file.path(out_dir, "coloc_eb_summary.tsv"), sep = "\t")
fwrite(rbindlist(timing_rows, fill = TRUE), file.path(out_dir, "timings.tsv"), sep = "\t")
fwrite(data.table(warning = unique(warnings_seen)), file.path(out_dir, "warnings.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

status <- readLines("/proc/self/status", warn = FALSE)
writeLines(grep("^(VmPeak|VmHWM|Threads):", status, value = TRUE), file.path(out_dir, "process_memory.txt"))
message("completed MAPK3 EB sensitivity: ", length(fit_files), " fits")
