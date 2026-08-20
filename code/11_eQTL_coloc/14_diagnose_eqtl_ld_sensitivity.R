#!/usr/bin/env Rscript

## reconcile every non-MAPK3 credible-set change in the raw-LD sensitivity run.

pilot_lib <- paste0(
    "/home/gpertea/work/R/spatialDLPFC_mdd_bpd/processed-data/11_eQTL_coloc/",
    "seurat/coloc/susie_exploratory/Rlib"
)
.libPaths(c(pilot_lib, .libPaths()))
suppressPackageStartupMessages({
    library(data.table)
    library(pgenlibr)
    library(susieR)
})
data.table::setDTthreads(1L)

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || length(args) > 2L) {
    stop("usage: 14_diagnose_eqtl_ld_sensitivity.R BASE [TARGET_ID]")
}
base <- normalizePath(args[[1L]], mustWork = TRUE)
manifest <- fread(file.path(base, "target_manifest.tsv"))
comparisons <- fread(file.path(base, "aggregate", "fit_sensitivity.tsv.gz"))
out_dir <- file.path(base, "diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_table <- file.path(out_dir, "eqtl_ld_sensitivity_reconciliation.tsv")
out_rds <- file.path(out_dir, "eqtl_ld_sensitivity_reconciliation.rds")
old_rows <- if (file.exists(out_table)) fread(out_table) else data.table()
old_fits <- if (file.exists(out_rds)) readRDS(out_rds) else list()
targets <- manifest[
    target_id %chin% comparisons[cs_identical == FALSE | max_abs_pip_diff > 0.1, target_id] &
        mapk3_pilot_complete == FALSE
]
if (length(args) == 2L) targets <- targets[target_id == args[[2L]]]

repo_root <- "/home/gpertea/work/R/spatialDLPFC_mdd_bpd"
processed <- file.path(repo_root, "processed-data", "11_eQTL_coloc", "seurat")
tqtl_in <- file.path(processed, "tqtl_in")
plink_prefix <- file.path(repo_root, "processed-data", "00_genotypes", "plink2", "merged_maf05")
pvar <- fread(file.path(base, "inputs", "target_pvar.tsv.gz"))

## load the union once to avoid reopening and rereading PGEN for each sensitivity case
objects <- setNames(lapply(targets$target_id, function(id) {
    readRDS(file.path(base, "results", paste0(id, ".rds")))
}), targets$target_id)
ids_by_target <- lapply(objects, function(x) x$variant_audit[is.na(exclusion_reason), variant_id])
union_ids <- unique(unlist(ids_by_target, use.names = FALSE))
pv <- pvar[variant_id %chin% union_ids]
setorder(pv, var_idx)
stopifnot(nrow(pv) == length(union_ids))

psam <- fread(paste0(plink_prefix, ".psam"))
sample_col <- intersect(c("#IID", "IID"), names(psam))[[1L]]
pvar_object <- NewPvar(paste0(plink_prefix, ".pvar"))
pgen <- NewPgen(paste0(plink_prefix, ".pgen"), pvar = pvar_object)
on.exit(try(ClosePgen(pgen), silent = TRUE), add = TRUE)
X_all <- ReadIntList(pgen, pv$var_idx)
rownames(X_all) <- as.character(psam[[sample_col]])
colnames(X_all) <- pv$variant_id

susie_args <- list(
    L = 10L, scaled_prior_variance = 0.2,
    estimate_prior_variance = TRUE, estimate_residual_variance = TRUE,
    z_method = "wald", coverage = 0.95, min_abs_corr = 0.5,
    max_iter = 1000L, refine = FALSE, verbose = FALSE, R_finite = FALSE
)
fit_rss <- function(z, R, n) {
    do.call(susie_rss, c(list(z = z, R = R, n = n), susie_args))
}
cs_ids <- function(fit, ids) {
    if (is.null(fit$sets$cs) || !length(fit$sets$cs)) return("")
    paste(sort(vapply(
        fit$sets$cs, function(i) paste(sort(ids[i]), collapse = ","), character(1L)
    )), collapse = ";")
}

rows <- list()
fit_objects <- list()
for (i in seq_len(nrow(targets))) {
    t <- targets[i]
    id <- t$target_id
    x <- objects[[id]]
    ids <- ids_by_target[[id]]
    cov_raw <- fread(file.path(tqtl_in, sprintf("%s.gene.covars.txt", t$dataset_id)))
    C <- t(as.matrix(cov_raw[, -1L]))
    storage.mode(C) <- "double"
    samples <- rownames(C)

    expr <- fread(file.path(tqtl_in, sprintf("%s.gene.expr.bed.gz", t$dataset_id)))
    yrow <- expr[ID == t$gene_id]
    expr_samples <- setdiff(names(expr), c("#Chr", "start", "end", "ID"))
    stopifnot(nrow(yrow) == 1L, identical(samples, expr_samples))
    y <- as.numeric(yrow[, ..expr_samples])

    X <- X_all[samples, ids, drop = FALSE]
    for (j in which(colSums(is.na(X)) > 0L)) X[is.na(X[, j]), j] <- mean(X[, j], na.rm = TRUE)
    qrc <- qr(C, tol = 1e-7, LAPACK = FALSE)
    Xr <- qr.resid(qrc, X)
    yr <- drop(qr.resid(qrc, y))
    R_rank <- cor(Xr)
    R_raw <- cor(X)
    dimnames(R_rank) <- dimnames(R_raw) <- list(ids, ids)

    ## recover residual-space degrees of freedom and exact compressed individual data
    n_effective <- length(samples) - qrc$rank
    xx <- colSums(Xr^2)
    xy <- drop(crossprod(Xr, yr))
    beta_rank <- xy / xx
    sse_rank <- sum(yr^2) - beta_rank^2 * xx
    se_rank <- sqrt((sse_rank / (n_effective - 1L)) / xx)
    z_rank <- beta_rank / se_rank
    names(z_rank) <- ids
    fit_rank <- fit_rss(z_rank, R_rank, n_effective)

    Q_complete <- qr.Q(qrc, complete = TRUE)
    U <- Q_complete[, seq.int(qrc$rank + 1L, nrow(C)), drop = FALSE]
    X_compressed <- crossprod(U, X)
    y_compressed <- drop(crossprod(U, y))
    fit_compressed <- susie(
        X_compressed, y_compressed, L = 10L, scaled_prior_variance = 0.2,
        estimate_prior_variance = TRUE, estimate_residual_variance = TRUE,
        intercept = FALSE, standardize = TRUE, coverage = 0.95,
        min_abs_corr = 0.5, max_iter = 1000L, refine = FALSE, verbose = FALSE
    )
    names(fit_compressed$pip) <- ids

    nominal <- fread(t$nominal_file)
    stored_z <- nominal[match(ids, variant_id), slope / slope_se]
    names(stored_z) <- ids
    rows[[id]] <- data.table(
        target_id = id, disorder = t$disorder, context = t$context,
        dataset_id = t$dataset_id, gene_id = t$gene_id, gene_name = t$gene_name,
        n = length(samples), covariate_rank = qrc$rank, n_effective = n_effective,
        n_variants = length(ids),
        max_abs_rank_vs_raw_LD = max(abs(R_rank - R_raw)),
        median_abs_rank_vs_raw_LD = median(abs(R_rank - R_raw)),
        cor_stored_rank_z = cor(stored_z, z_rank),
        max_abs_stored_rank_z = max(abs(stored_z - z_rank)),
        median_abs_stored_rank_z = median(abs(stored_z - z_rank)),
        max_abs_pip_primary_vs_raw = max(abs(x$fit_eqtl$pip - x$fit_eqtl_raw$pip)),
        max_abs_pip_primary_vs_rank_df = max(abs(x$fit_eqtl$pip - fit_rank$pip)),
        max_abs_pip_rank_df_vs_compressed = max(abs(fit_rank$pip - fit_compressed$pip)),
        n_cs_primary = length(x$fit_eqtl$sets$cs),
        n_cs_raw = length(x$fit_eqtl_raw$sets$cs),
        n_cs_rank_df = length(fit_rank$sets$cs),
        n_cs_compressed = length(fit_compressed$sets$cs),
        cs_primary = cs_ids(x$fit_eqtl, ids), cs_raw = cs_ids(x$fit_eqtl_raw, ids),
        cs_rank_df = cs_ids(fit_rank, ids), cs_compressed = cs_ids(fit_compressed, ids)
    )
    fit_objects[[id]] <- list(fit_rank_df = fit_rank, fit_compressed = fit_compressed)
    message(id, " complete")
}

all_rows <- rbindlist(c(list(old_rows), rows), fill = TRUE)
setorder(all_rows, target_id)
all_rows <- all_rows[!duplicated(target_id, fromLast = TRUE)]
fwrite(all_rows, out_table, sep = "\t")
saveRDS(c(old_fits, fit_objects), out_rds, compress = "gzip")
