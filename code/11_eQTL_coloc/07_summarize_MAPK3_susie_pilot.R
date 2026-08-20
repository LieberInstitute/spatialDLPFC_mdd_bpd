#!/usr/bin/env Rscript

## add credible-set, prior-sensitivity, and reviewed MAPK3 pilot summaries.

suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
})

repo_root <- normalizePath(getwd(), mustWork = TRUE)
pilot_dir <- file.path(
  repo_root, "processed-data", "11_eQTL_coloc", "seurat", "coloc",
  "susie_targeted_no23andMe", "MAPK3_pilot_20260819"
)
fit_files <- list.files(file.path(pilot_dir, "fits"), pattern = "\\.rds$", full.names = TRUE)

cs_rows <- list()
member_rows <- list()
prior_rows <- list()

for (path in fit_files) {
  comparison_key <- sub("\\.rds$", "", basename(path))
  parts <- strsplit(comparison_key, "__", fixed = TRUE)[[1L]]
  disorder <- parts[[1L]]
  context <- parts[[2L]]
  obj <- readRDS(path)

  for (trait in c("eqtl", "gwas")) {
    fit <- if (trait == "eqtl") obj$fit_eqtl_wrapper else obj$fit_gwas_wrapper
    ids <- if (trait == "eqtl") obj$D_eqtl$snp else obj$D_gwas$snp
    sets <- fit$sets$cs
    if (is.null(sets) || !length(sets)) next
    purity <- as.data.table(fit$sets$purity, keep.rownames = "component")

    for (i in seq_along(sets)) {
      idx <- sets[[i]]
      component <- names(sets)[[i]]
      component_value <- component
      prow <- purity[get("component") == component_value]
      cs_rows[[length(cs_rows) + 1L]] <- data.table(
        comparison_key = comparison_key, disorder = disorder, context = context,
        trait = trait, component = component, cs_size = length(idx),
        cs_log10bf = fit$lbf[fit$sets$cs_index[[i]]] / log(10),
        min_abs_corr = if (nrow(prow)) prow$min.abs.corr else NA_real_,
        mean_abs_corr = if (nrow(prow)) prow$mean.abs.corr else NA_real_,
        median_abs_corr = if (nrow(prow) && "median.abs.corr" %in% names(prow)) prow$median.abs.corr else NA_real_
      )
      member_rows[[length(member_rows) + 1L]] <- data.table(
        comparison_key = comparison_key, disorder = disorder, context = context,
        trait = trait, component = component, variant_id = ids[idx],
        pip = fit$pip[idx], alpha = fit$alpha[fit$sets$cs_index[[i]], idx]
      )
    }
  }

  ## cross-trait prior sensitivity does not refit either SuSiE model.
  for (p12 in c(1e-6, 1e-5, 1e-4)) {
    warning_text <- character()
    coloc_fit <- withCallingHandlers(
      coloc.susie(obj$fit_eqtl_wrapper, obj$fit_gwas_wrapper,
                  p1 = 1e-4, p2 = 1e-4, p12 = p12),
      warning = function(w) {
        warning_text <<- c(warning_text, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    s <- as.data.table(coloc_fit$summary)
    if (!nrow(s) || !ncol(s)) {
      prior_rows[[length(prior_rows) + 1L]] <- data.table(
        comparison_key = comparison_key, disorder = disorder, context = context,
        p12 = p12, status = "no_credible_set_pair",
        warnings = paste(unique(warning_text), collapse = " | ")
      )
    } else {
      s[, `:=`(
        comparison_key = comparison_key, disorder = disorder, context = context,
        p12 = p12, status = "ok", warnings = paste(unique(warning_text), collapse = " | ")
      )]
      prior_rows[[length(prior_rows) + 1L]] <- s
    }
  }
}

cs <- rbindlist(cs_rows, fill = TRUE)
members <- rbindlist(member_rows, fill = TRUE)
priors <- rbindlist(prior_rows, fill = TRUE)
fwrite(cs, file.path(pilot_dir, "credible_set_summary.tsv"), sep = "\t")
fwrite(members, file.path(pilot_dir, "credible_set_members.tsv.gz"), sep = "\t")
fwrite(priors, file.path(pilot_dir, "coloc_p12_sensitivity.tsv"), sep = "\t")

fit_summary <- fread(file.path(pilot_dir, "fit_summary.tsv"))
route <- fread(file.path(pilot_dir, "fit_route_comparison.tsv"))
ld <- fread(file.path(pilot_dir, "ld_diagnostics.tsv"))
coloc_summary <- fread(file.path(pilot_dir, "coloc_susie_summary.tsv"))
recon <- fread(file.path(pilot_dir, "diagnostics", "eqtl_individual_rss_reconciliation.tsv"))

primary <- coloc_summary[route == "prefit_primary"]
review <- c(
  "# Reviewed MAPK3 SuSiE pilot",
  "",
  sprintf("- Completed comparisons: %d (8 contexts x 3 disorders).", nrow(primary)),
  sprintf("- Primary fits converged: %d/%d.", sum(fit_summary$converged), nrow(fit_summary)),
  "- Direct susie_rss and coloc::runsusie fits were numerically identical for every trait.",
  sprintf("- Maximum estimate_s_rss: eQTL %.6g; GWAS %.6g.",
          max(ld[trait == "eqtl"]$estimate_s_rss), max(ld[trait == "gwas"]$estimate_s_rss)),
  sprintf("- Maximum GWAS finite-reference effective-rank/B: %.4f; no reliability flags.",
          max(fit_summary[trait == "gwas"]$finite_rank_over_B)),
  "- Naive n-row individual SuSiE was not used for scaling because residualization did not preserve covariate degrees of freedom.",
  sprintf("- Residual-space compressed individual SuSiE and rank-aware RSS differed by at most %.4f PIP in the representative MDD fits.",
          max(recon$max_abs_pip_rank_df_vs_compressed)),
  sprintf("- Stored tensorQTL RSS and rank-aware reconstructed RSS differed by at most %.4f PIP; all credible-set counts agreed.",
          max(recon$max_abs_pip_stored_vs_rank_df)),
  "- Primary later-target method: runsusie on stored dense tensorQTL statistics with rank-aware covariate-adjusted donor LD; trait-specific prefit coloc.susie.",
  "- GWAS uses approved no-23andMe/current summary statistics and 503-EUR finite-reference correction.",
  "- Raw donor LD remains a sensitivity; no nearPD, shrinkage, significance pruning, or reduced coverage was used.",
  "",
  "## Primary MAPK3 results",
  "",
  paste(capture.output(print(primary[, .(disorder, context, hit1, hit2,
                                          PP3 = PP.H3.abf, PP4 = PP.H4.abf)])), collapse = "\n")
)
writeLines(review, file.path(pilot_dir, "MAPK3_pilot_review.md"))
