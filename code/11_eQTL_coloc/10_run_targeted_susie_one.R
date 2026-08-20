#!/usr/bin/env Rscript

## run one checkpointed targeted SuSiE comparison from isolated approved inputs.
## coloc.abf is never called and LD matrices are not retained in checkpoints.

pilot_lib <- paste0(
    "/home/gpertea/work/R/spatialDLPFC_mdd_bpd/processed-data/11_eQTL_coloc/",
    "seurat/coloc/susie_exploratory/Rlib"
)
.libPaths(c(pilot_lib, .libPaths()))

suppressPackageStartupMessages({
    library(coloc)
    library(data.table)
    library(pgenlibr)
    library(susieR)
})

data.table::setDTthreads(1L)
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: 10_run_targeted_susie_one.R BASE TARGET_ID")
base <- normalizePath(args[[1L]], mustWork = TRUE)
target_key <- args[[2L]]

## use a private directory and expose one atomic checkpoint at completion
result_dir <- file.path(base, "results")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
final_file <- file.path(result_dir, paste0(target_key, ".rds"))
if (file.exists(final_file)) {
    message(target_key, " already complete")
    quit(save = "no", status = 0L)
}
tmp_file <- tempfile(paste0(".", target_key, "_"), tmpdir = result_dir, fileext = ".rds")
on.exit(unlink(tmp_file), add = TRUE)

manifest <- fread(file.path(base, "target_manifest.tsv"))
target <- manifest[target_id == target_key]
if (nrow(target) != 1L) stop("target ID not unique: ", target_key)
if (isTRUE(target$mapk3_pilot_complete)) stop("MAPK3 pilot target is handled separately")

repo_root <- "/home/gpertea/work/R/spatialDLPFC_mdd_bpd"
processed <- file.path(repo_root, "processed-data", "11_eQTL_coloc", "seurat")
tqtl_in <- file.path(processed, "tqtl_in")
plink_prefix <- file.path(repo_root, "processed-data", "00_genotypes", "plink2", "merged_maf05")
ref_sample_file <- file.path(
    processed, "coloc", "susie_exploratory", "reference",
    "1000G_phase3_EUR_2504-unrelated.samples"
)

timings <- list()
warnings_seen <- character()
timed <- function(step, expr) {
    gc()
    started <- Sys.time()
    tm <- system.time(value <- withCallingHandlers(
        force(expr),
        warning = function(w) {
            warnings_seen <<- c(warnings_seen, paste0(step, ": ", conditionMessage(w)))
            invokeRestart("muffleWarning")
        }
    ))
    timings[[length(timings) + 1L]] <<- data.table(
        step = step, started = format(started, "%Y-%m-%d %H:%M:%S %z"),
        elapsed_seconds = unname(tm[["elapsed"]]),
        user_seconds = unname(tm[["user.self"]]),
        system_seconds = unname(tm[["sys.self"]])
    )
    value
}

## load this gene's stored tensorQTL nominal statistics and exact PGEN rows
eqtl <- fread(target$nominal_file)
if (anyDuplicated(eqtl$variant_id)) stop("duplicate nominal variant IDs")
pvar <- fread(file.path(base, "inputs", "target_pvar.tsv.gz"))
pv <- pvar[variant_id %chin% eqtl$variant_id]
setorder(pv, var_idx)
if (nrow(pv) != nrow(eqtl)) stop("nominal/PVAR variant mismatch")

psam <- fread(paste0(plink_prefix, ".psam"))
sample_col <- intersect(c("#IID", "IID"), names(psam))[[1L]]
pvar_object <- NewPvar(paste0(plink_prefix, ".pvar"))
pgen <- NewPgen(paste0(plink_prefix, ".pgen"), pvar = pvar_object)
on.exit(try(ClosePgen(pgen), silent = TRUE), add = TRUE)
X_all <- timed("read_pgen_hardcalls", ReadIntList(pgen, pv$var_idx))
rownames(X_all) <- as.character(psam[[sample_col]])
colnames(X_all) <- pv$variant_id

## reproduce the documented expression/covariate sample order and rank-aware projection
cov_path <- file.path(tqtl_in, sprintf("%s.gene.covars.txt", target$dataset_id))
cov_raw <- fread(cov_path)
C <- t(as.matrix(cov_raw[, -1L]))
storage.mode(C) <- "double"
rownames(C) <- names(cov_raw)[-1L]
samples <- rownames(C)
if (!all(samples %chin% rownames(X_all))) stop("covariate/PGEN sample mismatch")
X <- X_all[samples, eqtl$variant_id, drop = FALSE]
missing_eqtl <- colSums(is.na(X))
means <- colMeans(X, na.rm = TRUE)
for (j in which(missing_eqtl > 0L)) X[is.na(X[, j]), j] <- means[[j]]
qrc <- qr(C, tol = 1e-07, LAPACK = FALSE)
X_resid <- timed("residualize_eqtl_genotypes", qr.resid(qrc, X))
rm(X_all)

## load fresh approved GWAS and exact 503-EUR reference subsets
gwas_manifest <- fread(file.path(base, "gwas_manifest.tsv"))
gm <- gwas_manifest[disorder == target$disorder]
gwas <- fread(gm$approved_file)
gwas <- gwas[variant_id %chin% eqtl$variant_id]
gwas_all <- fread(gm$harmonization_file)
gwas_all <- gwas_all[variant_id %chin% eqtl$variant_id]

ref_file <- file.path(
    base, "inputs", "reference",
    sprintf("%s.1000G_highcov.EUR_unrelated.GT.tsv.gz", target$chr)
)
ref <- fread(ref_file)
ref <- ref[variant_id %chin% eqtl$variant_id]
ref_samples <- setdiff(names(ref), "variant_id")
expected_ref_samples <- scan(ref_sample_file, what = "", quiet = TRUE)
if (!identical(ref_samples, expected_ref_samples)) stop("reference sample order mismatch")

## construct an explicit exclusion ledger before any LD or model fitting
audit <- eqtl[, .(
    variant_id, eqtl_beta = slope, eqtl_se = slope_se,
    eqtl_p = pval_nominal, eqtl_af = af
)]
audit <- merge(audit, gwas[, .(
    variant_id, rsid, gwas_beta = beta, gwas_se = beta_se,
    gwas_p = p, N, ncas, impinfo, match_mode
)], by = "variant_id", all.x = TRUE)
gwas_excl <- gwas_all[!is.na(exclusion_reason) & nzchar(exclusion_reason), .(
    variant_id, gwas_input_exclusion = exclusion_reason
)]
audit <- merge(audit, gwas_excl, by = "variant_id", all.x = TRUE)
audit[, exclusion_reason := NA_character_]
audit[!is.finite(eqtl_beta) | !is.finite(eqtl_se) | eqtl_se <= 0,
      exclusion_reason := "invalid_eqtl_stats"]
audit[is.na(exclusion_reason) & !is.na(gwas_input_exclusion),
      exclusion_reason := gwas_input_exclusion]
audit[is.na(exclusion_reason) & is.na(gwas_beta),
      exclusion_reason := "no_GWAS_allele_match"]
audit[is.na(exclusion_reason) & !variant_id %chin% ref$variant_id,
      exclusion_reason := "absent_1000G_reference"]

ids <- audit[is.na(exclusion_reason), variant_id]
Xg <- t(as.matrix(ref[match(ids, variant_id), ..ref_samples]))
storage.mode(Xg) <- "double"
colnames(Xg) <- ids
missing_ref <- colSums(is.na(Xg))
ref_means <- colMeans(Xg, na.rm = TRUE)
for (j in which(missing_ref > 0L)) Xg[is.na(Xg[, j]), j] <- ref_means[[j]]
bad_ref <- names(which(!is.finite(apply(Xg, 2L, sd)) | apply(Xg, 2L, sd) == 0))
audit[variant_id %chin% bad_ref, exclusion_reason := "monomorphic_1000G_reference"]

ids <- audit[is.na(exclusion_reason), variant_id]
bad_eqtl <- names(which(!is.finite(apply(X_resid[, ids, drop = FALSE], 2L, sd)) |
                        apply(X_resid[, ids, drop = FALSE], 2L, sd) == 0))
audit[variant_id %chin% bad_eqtl, exclusion_reason := "monomorphic_eqtl_residual"]
ids <- audit[is.na(exclusion_reason), variant_id]
if (length(ids) < 10L) {
    print(audit[, .N, by = .(exclusion_reason)][order(-N)])
    stop("too few variants after exclusions")
}

e <- eqtl[match(ids, variant_id)]
g <- gwas[match(ids, variant_id)]
Xe <- X_resid[, ids, drop = FALSE]
Xe_raw <- X[, ids, drop = FALSE]
Xg <- t(as.matrix(ref[match(ids, variant_id), ..ref_samples]))
storage.mode(Xg) <- "double"
colnames(Xg) <- ids
for (j in which(colSums(is.na(Xg)) > 0L)) Xg[is.na(Xg[, j]), j] <- mean(Xg[, j], na.rm = TRUE)
rm(ref, X, X_resid, pvar, pv, gwas_all)

## calculate exact low-rank LD compatibility diagnostics without a p-by-p eigen decomposition
lowrank_diagnostics <- function(z, Xs, n, trait) {
    p <- ncol(Xs)
    z_adj <- sqrt((n - 1) / (z^2 + n - 2)) * z
    K <- tcrossprod(Xs) / (nrow(Xs) - 1)
    ek <- eigen(K, symmetric = TRUE)
    keep <- which(ek$values > 1e-8)
    d <- ek$values[keep]
    U <- ek$vectors[, keep, drop = FALSE]
    V <- crossprod(Xs, U)
    V <- sweep(V, 2L, sqrt((nrow(Xs) - 1) * d), "/")
    proj <- drop(crossprod(V, z_adj))
    null_ss <- max(0, sum(z_adj^2) - sum(proj^2))
    objective <- function(s) {
        if (s <= 0) return(if (null_ss < 1e-10) 0.5 * sum(log(d) + proj^2 / d) else Inf)
        den <- (1 - s) * d + s
        0.5 * (sum(log(den) + proj^2 / den) + (p - length(d)) * log(s) + null_ss / s)
    }
    s <- optimize(objective, c(0, 1), tol = 1e-8)$minimum
    s_safe <- max(s, 1e-12)
    den <- (1 - s_safe) * d + s_safe
    delta <- 1 / den - 1 / s_safe
    pdiag <- 1 / s_safe + drop((V^2) %*% delta)
    pz <- z_adj / s_safe + drop(V %*% (delta * proj))
    z_std <- pz / sqrt(pdiag)
    data.table(
        trait = trait, n = n, n_variants = p, effective_rank = length(d),
        rank_over_variants = length(d) / p,
        min_eigenvalue = if (length(d) < p) 0 else min(d),
        max_eigenvalue = max(d), estimate_s_rss = s,
        kriging_max_abs_z_std_diff = max(abs(z_std)),
        kriging_n_abs_z_std_gt3 = sum(abs(z_std) > 3),
        diag_max_abs_error = max(abs(colSums(Xs^2) / (nrow(Xs) - 1) - 1))
    )
}

susie_common <- list(
    L = 10L, scaled_prior_variance = 0.2,
    estimate_prior_variance = TRUE, z_method = "wald",
    coverage = 0.95, min_abs_corr = 0.5,
    max_iter = 1000L, refine = FALSE, verbose = FALSE
)
susie_eqtl <- modifyList(susie_common, list(estimate_residual_variance = TRUE, R_finite = FALSE))
susie_gwas <- modifyList(susie_common, list(estimate_residual_variance = FALSE,
                                            R_finite = length(ref_samples)))
n_eqtl <- length(samples)
n_gwas <- median(g$N, na.rm = TRUE)
z_eqtl <- e$slope / e$slope_se
z_gwas <- g$beta / g$beta_se
names(z_eqtl) <- names(z_gwas) <- ids

## standardize once, fit the validated faster runsusie route, then release each LD matrix
Xe_std <- scale(Xe)
ld_eqtl <- timed("eqtl_ld_diagnostics", lowrank_diagnostics(z_eqtl, Xe_std, n_eqtl, "eqtl"))
R_eqtl <- timed("build_eqtl_ld", crossprod(Xe_std) / (nrow(Xe_std) - 1))
dimnames(R_eqtl) <- list(ids, ids)
D_eqtl <- list(
    beta = e$slope, varbeta = e$slope_se^2, z = z_eqtl, snp = ids,
    LD = R_eqtl, N = n_eqtl, MAF = pmin(e$af, 1 - e$af), type = "quant"
)
fit_eqtl <- timed("runsusie_eqtl_primary", do.call(coloc::runsusie, c(list(d = D_eqtl), susie_eqtl)))
rm(D_eqtl, R_eqtl)

Xe_raw_std <- scale(Xe_raw)
R_raw <- timed("build_eqtl_raw_ld", crossprod(Xe_raw_std) / (nrow(Xe_raw_std) - 1))
dimnames(R_raw) <- list(ids, ids)
D_raw <- list(
    beta = e$slope, varbeta = e$slope_se^2, z = z_eqtl, snp = ids,
    LD = R_raw, N = n_eqtl, MAF = pmin(e$af, 1 - e$af), type = "quant"
)
fit_eqtl_raw <- timed("runsusie_eqtl_raw_ld_sensitivity",
                      do.call(coloc::runsusie, c(list(d = D_raw), susie_eqtl)))
rm(D_raw, R_raw, Xe_raw_std, Xe_raw)

Xg_std <- scale(Xg)
ld_gwas <- timed("gwas_ld_diagnostics", lowrank_diagnostics(z_gwas, Xg_std, n_gwas, "gwas"))
R_gwas <- timed("build_gwas_ld", crossprod(Xg_std) / (nrow(Xg_std) - 1))
dimnames(R_gwas) <- list(ids, ids)
D_gwas <- list(
    beta = g$beta, varbeta = g$beta_se^2, z = z_gwas, snp = ids,
    LD = R_gwas, N = n_gwas, type = "cc",
    s = median(g$ncas / g$N, na.rm = TRUE)
)
fit_gwas <- timed("runsusie_gwas_primary", do.call(coloc::runsusie, c(list(d = D_gwas), susie_gwas)))
rm(D_gwas, R_gwas, Xg_std, Xg, Xe_std, Xe)

## coloc prefit objects retain trait-specific SuSiE settings; prior sensitivity is cheap
coloc_primary <- timed(
    "coloc_susie_primary",
    coloc::coloc.susie(fit_eqtl, fit_gwas, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
)
p12_values <- c(1e-6, 1e-5, 1e-4)
coloc_sensitivity <- lapply(p12_values, function(p12) {
    x <- withCallingHandlers(
        coloc::coloc.susie(fit_eqtl, fit_gwas, p1 = 1e-4, p2 = 1e-4, p12 = p12),
        warning = function(w) {
            warnings_seen <<- c(
                warnings_seen,
                paste0("coloc_susie_p12_", p12, ": ", conditionMessage(w))
            )
            invokeRestart("muffleWarning")
        }
    )
    s <- as.data.table(x$summary)
    if (!nrow(s)) s <- data.table(status = "no_credible_set_pair")
    s[, p12 := p12]
    s
})
coloc_sensitivity <- rbindlist(coloc_sensitivity, fill = TRUE)

canonical_cs <- function(fit) {
    cs <- fit$sets$cs
    if (is.null(cs) || !length(cs)) return(character())
    sort(vapply(cs, function(i) paste(sort(names(fit$pip)[i]), collapse = ","), character(1L)))
}
fit_comparison <- data.table(
    comparison = "eqtl_rank_adjusted_vs_raw_LD",
    max_abs_pip_diff = max(abs(fit_eqtl$pip - fit_eqtl_raw$pip)),
    cs_identical = identical(canonical_cs(fit_eqtl), canonical_cs(fit_eqtl_raw)),
    n_cs_primary = length(fit_eqtl$sets$cs), n_cs_sensitivity = length(fit_eqtl_raw$sets$cs)
)

summarize_fit <- function(fit, trait, n) {
    finite <- fit$R_finite_diagnostics
    data.table(
        trait = trait, n = n, n_variants = length(ids), converged = isTRUE(fit$converged),
        niter = fit$niter, n_credible_sets = length(fit$sets$cs), max_pip = max(fit$pip),
        top_variant = names(fit$pip)[which.max(fit$pip)], final_elbo = tail(fit$elbo, 1L),
        residual_variance = fit$sigma2,
        finite_effective_rank = if (is.null(finite)) NA_real_ else finite$effective_rank,
        finite_rank_over_B = if (is.null(finite)) NA_real_ else finite$r_over_B,
        finite_max_penalty = if (is.null(finite)) NA_real_ else max(finite$per_variable_penalty),
        finite_sensitivity_flag = if (is.null(finite)) NA else finite$R_sensitivity_flag,
        finite_reliability_flag = if (is.null(finite)) NA else finite$R_reliability_flag
    )
}
fit_summary <- rbindlist(list(
    summarize_fit(fit_eqtl, "eqtl", n_eqtl),
    summarize_fit(fit_gwas, "gwas", n_gwas)
))
exclusions <- audit[, .N, by = .(
    exclusion_reason = fifelse(is.na(exclusion_reason), "included", exclusion_reason)
)]

## read Linux high-water RSS for per-target memory accounting
status <- readLines("/proc/self/status", warn = FALSE)
vmhwm <- grep("^VmHWM:", status, value = TRUE)
peak_rss_kb <- as.numeric(sub("^VmHWM:[[:space:]]*([0-9]+).*$", "\\1", vmhwm))
metadata <- data.table(
    target_id = target_key, disorder = target$disorder, context = target$context,
    dataset_id = target$dataset_id, gene_id = target$gene_id, gene_name = target$gene_name,
    chr = target$chr, n_eqtl = n_eqtl, covariate_rank = qrc$rank,
    n_reference = length(ref_samples), n_variants = length(ids), peak_rss_kb = peak_rss_kb,
    completed = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
)

checkpoint <- list(
    metadata = metadata, variant_audit = audit, exclusions = exclusions,
    ld_diagnostics = rbindlist(list(ld_eqtl, ld_gwas), fill = TRUE),
    fit_summary = fit_summary, fit_comparison = fit_comparison,
    fit_eqtl = fit_eqtl, fit_eqtl_raw = fit_eqtl_raw, fit_gwas = fit_gwas,
    coloc_primary = coloc_primary, coloc_sensitivity = coloc_sensitivity,
    timings = rbindlist(timings), warnings = unique(warnings_seen),
    parameters = list(eqtl = susie_eqtl, gwas = susie_gwas, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
)
saveRDS(checkpoint, tmp_file, compress = "gzip")
if (!file.rename(tmp_file, final_file)) stop("failed to expose checkpoint")
message(target_key, " complete: ", length(ids), " variants")
