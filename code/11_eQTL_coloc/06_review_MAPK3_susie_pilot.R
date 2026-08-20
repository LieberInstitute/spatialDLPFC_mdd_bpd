#!/usr/bin/env Rscript

## reconcile individual-level and RSS eQTL fits for the MAPK3 pilot.

suppressPackageStartupMessages({
  library(data.table)
  library(pgenlibr)
  library(susieR)
})

repo_root <- normalizePath(file.path(getwd()), mustWork = TRUE)
if (!dir.exists(file.path(repo_root, ".git"))) stop("Run from repository root")
source(file.path(repo_root, "code", "11_eQTL_coloc", "utils.R"))

contexts <- c("astro", "inhb", "l2-3", "l4", "l5", "l6", "uvasc", "oligo")
gene_id <- "ENSG00000102882"
base_dir <- file.path(repo_root, "processed-data", "11_eQTL_coloc", "seurat")
pilot_dir <- file.path(base_dir, "coloc", "susie_targeted_no23andMe", "MAPK3_pilot_20260819")
tqtl_in <- file.path(base_dir, "tqtl_in")
plink_prefix <- file.path(repo_root, "processed-data", "00_genotypes", "plink2", "merged_maf05")

## use the MDD intersection as one representative variant set per context.
objects <- setNames(lapply(contexts, function(context) {
  readRDS(file.path(pilot_dir, "fits", paste0("MDD__", context, ".rds")))
}), contexts)
all_ids <- unique(unlist(lapply(objects, function(x) x$D_eqtl$snp), use.names = FALSE))

## recover PGEN indexes once and load the exact tensorQTL hardcalls.
pvar <- read_plink2_variant_table(plink_prefix)
pmap <- pvar[data.table(variant_id = all_ids), on = "variant_id", nomatch = 0L]
if (nrow(pmap) != length(all_ids)) stop("Missing PVAR variants")
setorder(pmap, var_idx)
psam <- fread(paste0(plink_prefix, ".psam"))
sample_col <- intersect(c("#IID", "IID"), names(psam))[[1L]]
pvar_object <- NewPvar(paste0(plink_prefix, ".pvar"))
pgen <- NewPgen(paste0(plink_prefix, ".pgen"), pvar = pvar_object)
X_all <- ReadIntList(pgen, pmap$var_idx)
rownames(X_all) <- psam[[sample_col]]
colnames(X_all) <- pmap$variant_id

args_rss <- list(
  L = 10L, scaled_prior_variance = 0.2,
  estimate_prior_variance = TRUE, estimate_residual_variance = TRUE,
  z_method = "wald", coverage = 0.95, min_abs_corr = 0.5,
  max_iter = 1000L, refine = FALSE, verbose = FALSE, R_finite = FALSE
)

fit_rss <- function(z, R, n) {
  do.call(susie_rss, c(list(z = z, R = R, n = n), args_rss))
}

cs_ids <- function(fit, ids) {
  if (is.null(fit$sets$cs) || !length(fit$sets$cs)) return("")
  paste(sort(vapply(fit$sets$cs, function(i) paste(sort(ids[i]), collapse = ","), character(1L))),
        collapse = ";")
}

rows <- list()
fits <- list()
for (context in contexts) {
  x <- objects[[context]]
  ids <- x$D_eqtl$snp

  cov_raw <- fread(file.path(tqtl_in, paste0(context, ".gene.covars.txt")))
  C <- t(as.matrix(cov_raw[, -1L]))
  storage.mode(C) <- "double"
  samples <- rownames(C)
  expr <- fread(file.path(tqtl_in, paste0(context, ".gene.expr.bed.gz")))
  yrow <- expr[ID == gene_id]
  expr_samples <- setdiff(names(expr), c("#Chr", "start", "end", "ID"))
  if (!identical(samples, expr_samples)) stop("Sample mismatch for ", context)
  y <- as.numeric(yrow[, ..expr_samples])

  X <- X_all[samples, ids, drop = FALSE]
  for (j in which(colSums(is.na(X)) > 0L)) X[is.na(X[, j]), j] <- mean(X[, j], na.rm = TRUE)
  qrc <- qr(C, tol = 1e-7, LAPACK = FALSE)
  Xr <- qr.resid(qrc, X)
  yr <- drop(qr.resid(qrc, y))
  Q_complete <- qr.Q(qrc, complete = TRUE)
  U <- Q_complete[, seq.int(qrc$rank + 1L, nrow(C)), drop = FALSE]
  X_compressed <- crossprod(U, X)
  y_compressed <- drop(crossprod(U, y))
  R <- cor(Xr)
  dimnames(R) <- list(ids, ids)

  ## reconstruct marginal statistics from the same rank-aware residual data.
  marginal <- univariate_regression(Xr, yr, center = FALSE, scale = FALSE)
  z_reconstructed <- marginal$betahat / marginal$sebetahat
  z_stored <- x$D_eqtl$z
  fit_reconstructed <- fit_rss(z_reconstructed, R, length(samples))
  fit_stored <- x$fit_eqtl_direct
  fit_individual <- x$fit_eqtl_individual

  ## form rank-aware marginal statistics with the correct residual degrees.
  n_effective <- length(samples) - qrc$rank
  fit_reconstructed_neff <- fit_rss(z_reconstructed, R, n_effective)
  xx <- colSums(Xr^2)
  xy <- drop(crossprod(Xr, yr))
  beta_rank <- xy / xx
  sse_rank <- sum(yr^2) - beta_rank^2 * xx
  df_rank <- n_effective - 1L
  se_rank <- sqrt((sse_rank / df_rank) / xx)
  z_rank_df <- beta_rank / se_rank
  fit_rank_df <- fit_rss(z_rank_df, R, n_effective)
  fit_compressed <- susie(
    X_compressed, y_compressed, L = 10L, scaled_prior_variance = 0.2,
    estimate_prior_variance = TRUE, estimate_residual_variance = TRUE,
    intercept = FALSE, standardize = TRUE, coverage = 0.95,
    min_abs_corr = 0.5, max_iter = 1000L, refine = FALSE, verbose = FALSE
  )
  names(fit_compressed$pip) <- ids
  fits[[context]] <- list(
    fit_stored = fit_stored, fit_reconstructed = fit_reconstructed,
    fit_reconstructed_neff = fit_reconstructed_neff,
    fit_rank_df = fit_rank_df, fit_compressed = fit_compressed,
    fit_individual = fit_individual, z_stored = z_stored,
    z_reconstructed = z_reconstructed, z_rank_df = z_rank_df,
    covariate_rank = qrc$rank
  )

  rows[[context]] <- data.table(
    context = context, n = length(samples), covariate_columns = ncol(C),
    covariate_rank = qrc$rank, n_effective = n_effective,
    max_abs_stored_vs_reconstructed_z = max(abs(z_stored - z_reconstructed)),
    median_abs_stored_vs_reconstructed_z = median(abs(z_stored - z_reconstructed)),
    cor_stored_reconstructed_z = cor(z_stored, z_reconstructed),
    max_abs_stored_vs_rank_df_z = max(abs(z_stored - z_rank_df)),
    median_abs_stored_vs_rank_df_z = median(abs(z_stored - z_rank_df)),
    cor_stored_rank_df_z = cor(z_stored, z_rank_df),
    max_pip_stored = max(fit_stored$pip),
    max_pip_reconstructed = max(fit_reconstructed$pip),
    max_pip_reconstructed_neff = max(fit_reconstructed_neff$pip),
    max_pip_rank_df = max(fit_rank_df$pip),
    max_pip_compressed = max(fit_compressed$pip),
    max_pip_individual = max(fit_individual$pip),
    max_abs_pip_stored_vs_reconstructed = max(abs(fit_stored$pip - fit_reconstructed$pip)),
    max_abs_pip_reconstructed_vs_individual = max(abs(fit_reconstructed$pip - fit_individual$pip)),
    max_abs_pip_reconstructed_neff_vs_individual = max(abs(fit_reconstructed_neff$pip - fit_individual$pip)),
    max_abs_pip_stored_vs_rank_df = max(abs(fit_stored$pip - fit_rank_df$pip)),
    max_abs_pip_rank_df_vs_compressed = max(abs(fit_rank_df$pip - fit_compressed$pip)),
    n_cs_stored = length(fit_stored$sets$cs),
    n_cs_reconstructed = length(fit_reconstructed$sets$cs),
    n_cs_reconstructed_neff = length(fit_reconstructed_neff$sets$cs),
    n_cs_rank_df = length(fit_rank_df$sets$cs),
    n_cs_compressed = length(fit_compressed$sets$cs),
    n_cs_individual = length(fit_individual$sets$cs),
    cs_stored = cs_ids(fit_stored, ids),
    cs_reconstructed = cs_ids(fit_reconstructed, ids),
    cs_reconstructed_neff = cs_ids(fit_reconstructed_neff, ids),
    cs_rank_df = cs_ids(fit_rank_df, ids),
    cs_compressed = cs_ids(fit_compressed, ids),
    cs_individual = cs_ids(fit_individual, ids)
  )
  message(context, " complete")
}

fwrite(rbindlist(rows, fill = TRUE),
       file.path(pilot_dir, "diagnostics", "eqtl_individual_rss_reconciliation.tsv"), sep = "\t")
saveRDS(fits, file.path(pilot_dir, "diagnostics", "eqtl_individual_rss_reconciliation.rds"),
        compress = "xz")
