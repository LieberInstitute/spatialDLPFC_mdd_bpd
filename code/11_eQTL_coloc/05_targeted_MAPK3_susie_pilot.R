#!/usr/bin/env Rscript

## run an isolated MAPK3 SuSiE/coloc pilot from approved source inputs.
## this script never reads production GWAS caches and never runs coloc.abf.

suppressPackageStartupMessages({
  library(arrow)
  library(coloc)
  library(data.table)
  library(pgenlibr)
  library(susieR)
})

data.table::setDTthreads(1L)
options(stringsAsFactors = FALSE)

## locate the repository without depending on the caller's working directory.
find_repo_root <- function(path = getwd()) {
  path <- normalizePath(path, mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(path, ".git"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Cannot locate repository root")
    path <- parent
  }
}

repo_root <- find_repo_root()
code_dir <- file.path(repo_root, "code", "11_eQTL_coloc")
source(file.path(code_dir, "utils.R"))

## environment filters support a small smoke run and the complete MAPK3 pilot.
split_env <- function(name, default) {
  value <- Sys.getenv(name, unset = paste(default, collapse = ","))
  trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
}

all_contexts <- c("astro", "inhb", "l2-3", "l4", "l5", "l6", "uvasc", "oligo")
all_disorders <- c("BD", "MDD", "SCZD")
contexts <- split_env("PILOT_CONTEXTS", all_contexts)
disorders <- split_env("PILOT_DISORDERS", all_disorders)
stopifnot(all(contexts %in% all_contexts), all(disorders %in% all_disorders))

label <- Sys.getenv("PILOT_LABEL", unset = "MAPK3_pilot_20260819")
processed_dir <- file.path(repo_root, "processed-data", "11_eQTL_coloc", "seurat")
out_dir <- file.path(processed_dir, "coloc", "susie_targeted_no23andMe", label)
input_dir <- file.path(out_dir, "inputs")
fit_dir <- file.path(out_dir, "fits")
diag_dir <- file.path(out_dir, "diagnostics")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

tqtl_in_dir <- file.path(processed_dir, "tqtl_in")
tqtl_out_dir <- file.path(processed_dir, "tqtl_out")
plink_prefix <- file.path(repo_root, "processed-data", "00_genotypes", "plink2", "merged_maf05")
ref_dir <- file.path(processed_dir, "coloc", "susie_exploratory", "reference")
ref_gt_file <- file.path(ref_dir, "MAPK3.1000G_highcov.EUR_unrelated.GT.tsv.gz")
ref_sample_file <- file.path(ref_dir, "1000G_phase3_EUR_2504-unrelated.samples")
tensorqtl_q_dir <- Sys.getenv("PILOT_TENSORQTL_Q_DIR", unset = "")

gene_id <- "ENSG00000102882"
gene_name <- "MAPK3"
si_min <- 0.8

## preserve per-step elapsed time in a simple append-only table.
timing_rows <- list()
timed <- function(step, key, expr) {
  gc()
  started <- Sys.time()
  tm <- system.time(value <- force(expr))
  timing_rows[[length(timing_rows) + 1L]] <<- data.table(
    step = step, comparison_key = key,
    started = format(started, "%Y-%m-%d %H:%M:%S %z"),
    elapsed_seconds = unname(tm[["elapsed"]]), user_seconds = unname(tm[["user.self"]]),
    system_seconds = unname(tm[["sys.self"]])
  )
  value
}

## return a stable SHA-256 without reading file contents into R.
sha256_file <- function(path) {
  out <- system2("sha256sum", shQuote(path), stdout = TRUE)
  sub("[[:space:]].*$", "", out[[1L]])
}

## read the approved all-donor manifest and verify requested sample counts.
manifest <- fread(file.path(tqtl_in_dir, "prep_manifest.csv"))
manifest <- manifest[dataset_id %in% contexts & split == "all"]
if (nrow(manifest) != length(contexts)) stop("Incomplete all-donor manifest")
manifest[, context_order := match(dataset_id, contexts)]
setorder(manifest, context_order)

## read dense MAPK3 nominal statistics from the chromosome parquet files.
read_eqtl <- function(context) {
  path <- file.path(tqtl_out_dir, sprintf("%s.gene.cis_qtl_pairs.chr16.parquet", context))
  cols <- c("phenotype_id", "variant_id", "start_distance", "af", "ma_samples",
            "ma_count", "pval_nominal", "slope", "slope_se")
  x <- as.data.table(read_parquet(path, col_select = tidyselect::all_of(cols)))
  x <- x[phenotype_id == gene_id]
  if (!nrow(x)) stop("MAPK3 missing from ", path)
  x[, pos := as.integer(sub("^[^:]+:([0-9]+):.*$", "\\1", variant_id))]
  setorder(x, pos, variant_id)
  if (anyDuplicated(x$variant_id)) stop("Duplicate eQTL variants for ", context)
  x
}

eqtl_by_context <- setNames(vector("list", length(contexts)), contexts)
for (context in contexts) {
  eqtl_by_context[[context]] <- timed("read_nominal_parquet", context, read_eqtl(context))
}
union_ids <- unique(unlist(lapply(eqtl_by_context, `[[`, "variant_id"), use.names = FALSE))

## load the PVAR once to recover one-based PGEN indexes and allele truth.
pvar <- timed("read_pvar", "MAPK3", read_plink2_variant_table(plink_prefix))
pvar_map <- pvar[data.table(variant_id = union_ids), on = "variant_id", nomatch = 0L]
if (nrow(pvar_map) != length(union_ids)) {
  stop(length(union_ids) - nrow(pvar_map), " nominal variants are absent from PVAR")
}
setorder(pvar_map, var_idx)

## load hardcalls exactly as the original tensorQTL PGEN reader did.
psam <- fread(paste0(plink_prefix, ".psam"))
sample_col <- intersect(c("#IID", "IID"), names(psam))[[1L]]
pvar_object <- NewPvar(paste0(plink_prefix, ".pvar"))
pgen <- NewPgen(paste0(plink_prefix, ".pgen"), pvar = pvar_object)
X_all <- timed("read_pgen_hardcalls", "MAPK3", ReadIntList(pgen, pvar_map$var_idx))
rownames(X_all) <- as.character(psam[[sample_col]])
colnames(X_all) <- pvar_map$variant_id

## mean-impute within each context after selecting its exact donor set.
read_context_data <- function(context) {
  cov_path <- file.path(tqtl_in_dir, sprintf("%s.gene.covars.txt", context))
  expr_path <- file.path(tqtl_in_dir, sprintf("%s.gene.expr.bed.gz", context))
  cov_raw <- fread(cov_path)
  cov_names <- as.character(cov_raw[[1L]])
  C <- t(as.matrix(cov_raw[, -1L]))
  storage.mode(C) <- "double"
  colnames(C) <- cov_names
  samples <- rownames(C)

  expr <- fread(expr_path)
  yrow <- expr[ID == gene_id]
  if (nrow(yrow) != 1L) stop("Expected one MAPK3 expression row for ", context)
  expr_samples <- setdiff(names(expr), c("#Chr", "start", "end", "ID"))
  if (!identical(samples, expr_samples)) stop("Expression/covariate order mismatch for ", context)
  if (!all(samples %in% rownames(X_all))) stop("PGEN sample mismatch for ", context)
  y <- as.numeric(yrow[, ..expr_samples])
  names(y) <- expr_samples

  ids <- eqtl_by_context[[context]]$variant_id
  X <- X_all[samples, ids, drop = FALSE]
  missing_before <- colSums(is.na(X))
  means <- colMeans(X, na.rm = TRUE)
  for (j in which(missing_before > 0L)) X[is.na(X[, j]), j] <- means[[j]]

  ## qr.resid retains the intended intercept and drops redundant columns by rank.
  qrc <- qr(C, tol = 1e-07, LAPACK = FALSE)
  X_resid <- qr.resid(qrc, X)
  y_resid <- drop(qr.resid(qrc, y))
  R_rank <- cor(X_resid)
  dimnames(R_rank) <- list(ids, ids)

  ## optional CPU PyTorch QR reproduces tensorQTL's rank-deficient projection.
  R_tensorqtl <- NULL
  q_path <- file.path(tensorqtl_q_dir, sprintf("%s.tensorqtl_Q.tsv.gz", context))
  if (nzchar(tensorqtl_q_dir) && file.exists(q_path)) {
    q_dt <- fread(q_path)
    q_samples <- as.character(q_dt[[1L]])
    Q <- as.matrix(q_dt[, -1L])
    storage.mode(Q) <- "double"
    if (!identical(q_samples, samples)) stop("tensorQTL Q sample mismatch for ", context)
    X_center <- scale(X, center = TRUE, scale = FALSE)
    X_tensorqtl <- X_center - Q %*% crossprod(Q, X_center)
    R_tensorqtl <- cor(X_tensorqtl)
    dimnames(R_tensorqtl) <- list(ids, ids)
  }

  list(
    samples = samples, covariates = C, covariate_rank = qrc$rank,
    X = X, X_resid = X_resid, y = y, y_resid = y_resid,
    R_rank = R_rank, R_tensorqtl = R_tensorqtl, missing_before = missing_before
  )
}

context_data <- setNames(vector("list", length(contexts)), contexts)
for (context in contexts) {
  context_data[[context]] <- timed("build_eqtl_rank_adjusted_ld", context, read_context_data(context))
}

## query each approved BCF freshly, then harmonize to PVAR REF/ALT.
query_gwas <- function(disorder) {
  bcf <- gwas_bcf_path(disorder, repo_root = repo_root)
  region_file <- file.path(input_dir, sprintf("MAPK3_%s_positions.bed", disorder))
  query_file <- tempfile(sprintf("MAPK3-%s-", disorder), fileext = ".tsv")
  err_file <- tempfile(sprintf("MAPK3-%s-", disorder), fileext = ".log")
  on.exit(unlink(c(query_file, err_file)), add = TRUE)

  regions <- unique(pvar_map[, .(CHROM, start0 = POS - 1L, end1 = POS)])
  fwrite(regions, region_file, sep = "\t", col.names = FALSE)
  status <- system2(
    "bcftools",
    c("query", "-R", shQuote(region_file),
      "-f", shQuote("%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%ES\t%SE\t%LP\t%NE\t%NS\t%NC\t%SI]\n"),
      shQuote(bcf)),
    stdout = query_file, stderr = err_file
  )
  if (status != 0L) stop(paste(readLines(err_file, warn = FALSE), collapse = "\n"))
  raw <- if (file.info(query_file)$size > 0) fread(query_file, header = FALSE) else data.table()
  if (!nrow(raw)) stop("No MAPK3 rows returned for ", disorder)
  setnames(raw, c("chr", "pos", "rsid", "a0", "a1", "beta", "beta_se",
                  "lp", "N", "ns", "ncas", "impinfo"))
  raw[, `:=`(
    chr = as.character(chr), pos = as.integer(pos),
    a0 = as.character(a0), a1 = as.character(a1),
    beta = suppressWarnings(as.numeric(beta)), beta_se = suppressWarnings(as.numeric(beta_se)),
    lp = suppressWarnings(as.numeric(lp)), N = suppressWarnings(as.numeric(N)),
    ncas = suppressWarnings(as.numeric(ncas)), impinfo = suppressWarnings(as.numeric(impinfo))
  )]
  raw[, p := 10^(-lp)]
  raw[, raw_row := .I]

  exact <- merge(raw, pvar_map, by.x = c("chr", "pos", "a0", "a1"),
                 by.y = c("CHROM", "POS", "REF", "ALT"), allow.cartesian = TRUE)
  exact[, `:=`(match_mode = "exact", match_rank = 1L)]
  swapped <- merge(raw, pvar_map, by.x = c("chr", "pos", "a0", "a1"),
                   by.y = c("CHROM", "POS", "ALT", "REF"), allow.cartesian = TRUE)
  swapped[, `:=`(beta = -beta, match_mode = "swapped", match_rank = 2L)]
  matched <- rbindlist(list(exact, swapped), use.names = TRUE, fill = TRUE)
  setorder(matched, variant_id, match_rank, p)
  matched <- matched[!duplicated(variant_id)]
  matched[, exclusion_reason := fcase(
    !is.finite(impinfo), "missing_SI",
    impinfo < si_min, "SI_below_0.8",
    !is.finite(beta) | !is.finite(beta_se) | beta_se <= 0, "invalid_beta_or_se",
    !is.finite(N) | N <= 0, "invalid_N",
    !is.finite(p) | p <= 0 | p > 1, "invalid_p",
    default = NA_character_
  )]
  keep <- matched[is.na(exclusion_reason)]
  fwrite(matched, file.path(input_dir, sprintf("MAPK3_%s_GWAS_harmonization_all.tsv.gz", disorder)), sep = "\t")
  fwrite(keep, file.path(input_dir, sprintf("MAPK3_%s_GWAS_approved.tsv.gz", disorder)), sep = "\t")
  list(data = keep, all = matched, bcf = bcf, raw_n = nrow(raw))
}

gwas_by_disorder <- setNames(vector("list", length(disorders)), disorders)
for (disorder in disorders) {
  gwas_by_disorder[[disorder]] <- timed("extract_and_harmonize_gwas", disorder, query_gwas(disorder))
}

## validate and decode the candidate 503-sample high-coverage reference table.
ref <- timed("read_reference_gt", "1000G_EUR", fread(ref_gt_file))
ref_ids <- as.character(ref$variant_id)
ref_samples <- setdiff(names(ref), "variant_id")
expected_ref_samples <- scan(ref_sample_file, what = character(), quiet = TRUE)
if (!setequal(ref_samples, expected_ref_samples) || anyDuplicated(ref_samples)) {
  stop("1000 Genomes reference sample list mismatch")
}

gt_chr <- as.matrix(ref[, ..ref_samples])
gt_lookup <- c("0|0" = 0, "0/0" = 0, "0|1" = 1, "1|0" = 1,
               "0/1" = 1, "1/0" = 1, "1|1" = 2, "1/1" = 2,
               ".|." = NA, "./." = NA, "." = NA)
gt_ref <- matrix(unname(gt_lookup[gt_chr]), nrow = nrow(gt_chr), ncol = ncol(gt_chr),
                 dimnames = list(ref_ids, ref_samples))
if (any(is.na(gt_ref) & !gt_chr %in% names(gt_lookup))) stop("Unexpected reference GT encoding")

## produce an auditable final variant set separately for every comparison.
build_pair <- function(context, disorder) {
  eqtl <- copy(eqtl_by_context[[context]])
  gwas <- gwas_by_disorder[[disorder]]$data
  gwas_all <- gwas_by_disorder[[disorder]]$all
  audit <- eqtl[, .(variant_id, pos, eqtl_beta = slope, eqtl_se = slope_se,
                    eqtl_p = pval_nominal, eqtl_af = af)]
  audit <- merge(audit, gwas[, .(variant_id, rsid, gwas_beta = beta,
                                 gwas_se = beta_se, gwas_p = p, N, ncas,
                                 impinfo, match_mode)], by = "variant_id", all.x = TRUE)
  gwas_exclusions <- gwas_all[!is.na(exclusion_reason), .(
    variant_id, gwas_input_exclusion = exclusion_reason
  )]
  audit <- merge(audit, gwas_exclusions, by = "variant_id", all.x = TRUE)
  audit[, exclusion_reason := NA_character_]
  audit[!is.finite(eqtl_beta) | !is.finite(eqtl_se) | eqtl_se <= 0,
        exclusion_reason := "invalid_eqtl_stats"]
  audit[is.na(exclusion_reason) & !is.na(gwas_input_exclusion),
        exclusion_reason := gwas_input_exclusion]
  audit[is.na(exclusion_reason) & is.na(gwas_beta),
        exclusion_reason := "no_GWAS_allele_match"]
  audit[is.na(exclusion_reason) & !variant_id %in% ref_ids,
        exclusion_reason := "absent_1000G_reference"]

  ids <- audit[is.na(exclusion_reason), variant_id]
  Xg <- t(gt_ref[ids, , drop = FALSE])
  missing_ref <- colSums(is.na(Xg))
  means <- colMeans(Xg, na.rm = TRUE)
  for (j in which(missing_ref > 0L)) Xg[is.na(Xg[, j]), j] <- means[[j]]
  ref_sd <- apply(Xg, 2L, sd)
  bad_ref <- names(ref_sd)[!is.finite(ref_sd) | ref_sd == 0]
  audit[variant_id %in% bad_ref, exclusion_reason := "monomorphic_1000G_reference"]

  ids <- audit[is.na(exclusion_reason), variant_id]
  Xe <- context_data[[context]]$X_resid[, ids, drop = FALSE]
  eqtl_sd <- apply(Xe, 2L, sd)
  bad_eqtl <- names(eqtl_sd)[!is.finite(eqtl_sd) | eqtl_sd == 0]
  audit[variant_id %in% bad_eqtl, exclusion_reason := "monomorphic_eqtl_residual"]
  ids <- audit[is.na(exclusion_reason), variant_id]
  if (length(ids) < 10L) stop("Too few variants for ", context, " / ", disorder)

  e <- eqtl[match(ids, variant_id)]
  g <- gwas[match(ids, variant_id)]
  Xe_raw <- context_data[[context]]$X[, ids, drop = FALSE]
  Xe_resid <- context_data[[context]]$X_resid[, ids, drop = FALSE]
  Xg <- t(gt_ref[ids, , drop = FALSE])
  for (j in which(colSums(is.na(Xg)) > 0L)) Xg[is.na(Xg[, j]), j] <- mean(Xg[, j], na.rm = TRUE)
  R_eqtl <- cor(Xe_resid)
  R_eqtl_raw <- cor(Xe_raw)
  R_eqtl_tensorqtl <- context_data[[context]]$R_tensorqtl
  if (!is.null(R_eqtl_tensorqtl)) R_eqtl_tensorqtl <- R_eqtl_tensorqtl[ids, ids, drop = FALSE]
  R_gwas <- cor(Xg)
  dimnames(R_eqtl) <- dimnames(R_eqtl_raw) <- dimnames(R_gwas) <- list(ids, ids)

  list(ids = ids, audit = audit, eqtl = e, gwas = g,
       Xe_resid = Xe_resid, y_resid = context_data[[context]]$y_resid,
       R_eqtl = R_eqtl, R_eqtl_raw = R_eqtl_raw,
       R_eqtl_tensorqtl = R_eqtl_tensorqtl, R_gwas = R_gwas,
       ref_missing = colSums(is.na(t(gt_ref[ids, , drop = FALSE]))))
}

## common arguments are used only for the dataset-route implementation check.
susie_args_common <- list(
  L = 10L, scaled_prior_variance = 0.2,
  estimate_prior_variance = TRUE, estimate_residual_variance = FALSE,
  z_method = "wald", coverage = 0.95, min_abs_corr = 0.5,
  max_iter = 1000L, refine = FALSE, verbose = FALSE
)
susie_args_eqtl <- modifyList(
  susie_args_common,
  list(estimate_residual_variance = TRUE, R_finite = FALSE)
)
susie_args_gwas <- modifyList(
  susie_args_common,
  list(R_finite = length(ref_samples))
)

do_susie_rss <- function(z, R, n, args) {
  do.call(susieR::susie_rss, c(list(z = z, R = R, n = n), args))
}

do_runsusie <- function(dataset, args) {
  do.call(coloc::runsusie, c(list(d = dataset), args))
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

## summarize signed-LD and summary-statistic compatibility diagnostics.
ld_diagnostics <- function(z, R, n, trait, key) {
  eigen_full <- eigen(R, symmetric = TRUE)
  eig <- eigen_full$values
  attr(R, "eigen") <- eigen_full
  s <- estimate_s_rss(z, R, n)
  kr <- kriging_rss(z, R, n, s = s)$conditional_dist
  data.table(
    comparison_key = key, trait = trait, n = n, n_variants = length(z),
    symmetry_max_abs = max(abs(R - t(R))), diag_max_abs_error = max(abs(diag(R) - 1)),
    min_eigenvalue = min(eig), n_eigen_below_minus_1e8 = sum(eig < -1e-8),
    estimate_s_rss = s, kriging_max_abs_z_std_diff = max(abs(kr$z_std_diff)),
    kriging_n_logLR_gt2_absz_gt2 = sum(kr$logLR > 2 & abs(kr$z) > 2)
  )
}

fit_rows <- list()
fit_compare_rows <- list()
coloc_rows <- list()
coloc_compare_rows <- list()
ld_rows <- list()
exclusion_rows <- list()

for (context in contexts) {
  for (disorder in disorders) {
    key <- paste(disorder, context, sep = "__")
    message(format(Sys.time()), " | starting ", key)
    pair <- timed("build_pair_and_ld", key, build_pair(context, disorder))
    fwrite(pair$audit, file.path(diag_dir, paste0(key, "_variant_audit.tsv.gz")), sep = "\t")
    exclusion_count <- pair$audit[, .N, by = .(
      exclusion_reason = fifelse(is.na(exclusion_reason), "included", exclusion_reason)
    )]
    exclusion_count[, `:=`(disorder = disorder, context = context)]
    setcolorder(exclusion_count, c("disorder", "context", "exclusion_reason", "N"))
    exclusion_rows[[key]] <- exclusion_count

    n_eqtl <- length(context_data[[context]]$samples)
    n_gwas <- median(pair$gwas$N, na.rm = TRUE)
    z_eqtl <- pair$eqtl$slope / pair$eqtl$slope_se
    z_gwas <- pair$gwas$beta / pair$gwas$beta_se
    names(z_eqtl) <- names(z_gwas) <- pair$ids
    D_eqtl <- list(
      beta = pair$eqtl$slope, varbeta = pair$eqtl$slope_se^2,
      z = z_eqtl, snp = pair$ids, LD = pair$R_eqtl, N = n_eqtl,
      MAF = pmin(pair$eqtl$af, 1 - pair$eqtl$af), type = "quant"
    )
    D_gwas <- list(
      beta = pair$gwas$beta, varbeta = pair$gwas$beta_se^2,
      z = z_gwas, snp = pair$ids, LD = pair$R_gwas, N = n_gwas, type = "cc",
      s = median(pair$gwas$ncas / pair$gwas$N, na.rm = TRUE)
    )

    ld_rows[[paste0(key, "_eqtl")]] <- timed(
      "ld_diagnostics", paste0(key, "__eqtl"),
      ld_diagnostics(z_eqtl, pair$R_eqtl, n_eqtl, "eqtl", key)
    )
    ld_rows[[paste0(key, "_gwas")]] <- timed(
      "ld_diagnostics", paste0(key, "__gwas"),
      ld_diagnostics(z_gwas, pair$R_gwas, n_gwas, "gwas", key)
    )

    fit_eqtl_direct <- timed(
      "susie_rss_direct_primary", paste0(key, "__eqtl"),
      do_susie_rss(z_eqtl, pair$R_eqtl, n_eqtl, susie_args_eqtl)
    )
    fit_eqtl_wrapper <- timed(
      "runsusie_primary", paste0(key, "__eqtl"),
      do_runsusie(D_eqtl, susie_args_eqtl)
    )
    fit_gwas_direct <- timed(
      "susie_rss_direct_primary", paste0(key, "__gwas"),
      do_susie_rss(z_gwas, pair$R_gwas, n_gwas, susie_args_gwas)
    )
    fit_gwas_wrapper <- timed(
      "runsusie_primary", paste0(key, "__gwas"),
      do_runsusie(D_gwas, susie_args_gwas)
    )

    ## individual-data SuSiE tests whether donor-level and RSS fits are reconcilable.
    fit_eqtl_individual <- timed(
      "susie_individual", paste0(key, "__eqtl"),
      susieR::susie(
        pair$Xe_resid, pair$y_resid, L = 10L, scaled_prior_variance = 0.2,
        estimate_prior_variance = TRUE, estimate_residual_variance = TRUE,
        intercept = FALSE, standardize = TRUE, coverage = 0.95,
        min_abs_corr = 0.5, max_iter = 1000L, refine = FALSE, verbose = FALSE
      )
    )
    names(fit_eqtl_individual$pip) <- pair$ids

    fit_eqtl_raw_ld <- timed(
      "susie_rss_raw_eqtl_ld_sensitivity", paste0(key, "__eqtl"),
      do_susie_rss(z_eqtl, pair$R_eqtl_raw, n_eqtl, susie_args_eqtl)
    )
    fit_eqtl_tensorqtl_ld <- NULL
    if (!is.null(pair$R_eqtl_tensorqtl)) {
      fit_eqtl_tensorqtl_ld <- timed(
        "susie_rss_tensorqtl_q_ld_sensitivity", paste0(key, "__eqtl"),
        do_susie_rss(z_eqtl, pair$R_eqtl_tensorqtl, n_eqtl, susie_args_eqtl)
      )
    }

    fit_compare_rows[[paste0(key, "_eqtl_wrapper")]] <- cbind(
      data.table(comparison_key = key, comparison = "eqtl_susie_rss_vs_runsusie"),
      compare_fits(fit_eqtl_direct, fit_eqtl_wrapper)
    )
    fit_compare_rows[[paste0(key, "_gwas_wrapper")]] <- cbind(
      data.table(comparison_key = key, comparison = "gwas_susie_rss_vs_runsusie"),
      compare_fits(fit_gwas_direct, fit_gwas_wrapper)
    )
    fit_compare_rows[[paste0(key, "_eqtl_individual")]] <- cbind(
      data.table(comparison_key = key, comparison = "eqtl_susie_rss_vs_individual_susie"),
      compare_fits(fit_eqtl_direct, fit_eqtl_individual)
    )
    fit_compare_rows[[paste0(key, "_eqtl_raw_ld")]] <- cbind(
      data.table(comparison_key = key, comparison = "eqtl_rank_adjusted_vs_raw_LD"),
      compare_fits(fit_eqtl_direct, fit_eqtl_raw_ld)
    )
    if (!is.null(fit_eqtl_tensorqtl_ld)) {
      fit_compare_rows[[paste0(key, "_eqtl_tensorqtl_ld")]] <- cbind(
        data.table(comparison_key = key, comparison = "eqtl_rank_adjusted_vs_tensorQTL_Q_LD"),
        compare_fits(fit_eqtl_direct, fit_eqtl_tensorqtl_ld)
      )
    }

    ## primary coloc uses trait-appropriate prefit objects.
    coloc_prefit_primary <- timed(
      "coloc_susie_prefit_primary", key,
      coloc::coloc.susie(fit_eqtl_wrapper, fit_gwas_wrapper,
                         p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
    )

    ## coloc's dataset route accepts only one common argument list for both traits.
    fit_eqtl_common <- timed(
      "runsusie_common_route_check", paste0(key, "__eqtl"),
      do_runsusie(D_eqtl, susie_args_common)
    )
    fit_gwas_common <- timed(
      "runsusie_common_route_check", paste0(key, "__gwas"),
      do_runsusie(D_gwas, susie_args_common)
    )
    coloc_prefit_common <- timed(
      "coloc_susie_prefit_common_route_check", key,
      coloc::coloc.susie(fit_eqtl_common, fit_gwas_common,
                         p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
    )
    coloc_dataset <- timed(
      "coloc_susie_dataset_common_route_check", key,
      coloc::coloc.susie(D_eqtl, D_gwas, susie.args = susie_args_common,
                         p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
    )

    summary_table <- function(x, route) {
      s <- as.data.table(x$summary)
      if (!ncol(s) || !nrow(s)) s <- data.table(status = "no_credible_set_pair")
      s[, `:=`(route = route, comparison_key = key,
               disorder = disorder, context = context)]
      s
    }
    s_primary <- summary_table(coloc_prefit_primary, "prefit_primary")
    s1 <- summary_table(coloc_prefit_common, "prefit_common_route_check")
    s2 <- summary_table(coloc_dataset, "dataset_common_route_check")
    coloc_rows[[key]] <- rbindlist(list(s_primary, s1, s2), use.names = TRUE, fill = TRUE)

    route_equal <- isTRUE(all.equal(coloc_prefit_common$summary, coloc_dataset$summary,
                                    tolerance = 1e-12, check.attributes = FALSE))
    pp_cols <- intersect(
      c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf"),
      intersect(names(s1), names(s2))
    )
    diffs <- unlist(lapply(pp_cols, function(x) abs(s1[[x]] - s2[[x]])))
    coloc_compare_rows[[key]] <- data.table(
      comparison_key = key, n_rows_prefit = nrow(s1), n_rows_dataset = nrow(s2),
      route_outputs_identical = route_equal,
      max_abs_pp_diff = if (length(diffs)) max(diffs, na.rm = TRUE) else NA_real_
    )

    saveRDS(
      list(
        ids = pair$ids, D_eqtl = D_eqtl, D_gwas = D_gwas,
        fit_eqtl_direct = fit_eqtl_direct, fit_eqtl_wrapper = fit_eqtl_wrapper,
        fit_eqtl_individual = fit_eqtl_individual,
        fit_eqtl_raw_ld = fit_eqtl_raw_ld,
        fit_eqtl_tensorqtl_ld = fit_eqtl_tensorqtl_ld,
        fit_gwas_direct = fit_gwas_direct, fit_gwas_wrapper = fit_gwas_wrapper,
        fit_eqtl_common = fit_eqtl_common, fit_gwas_common = fit_gwas_common,
        coloc_prefit_primary = coloc_prefit_primary,
        coloc_prefit_common = coloc_prefit_common, coloc_dataset_common = coloc_dataset,
        susie_args_eqtl = susie_args_eqtl, susie_args_gwas = susie_args_gwas,
        susie_args_common = susie_args_common
      ),
      file.path(fit_dir, paste0(key, ".rds")), compress = "xz"
    )

    for (trait in c("eqtl", "gwas")) {
      f <- if (trait == "eqtl") fit_eqtl_wrapper else fit_gwas_wrapper
      finite <- f$R_finite_diagnostics
      fit_rows[[paste(key, trait)]] <- data.table(
        comparison_key = key, disorder = disorder, context = context, trait = trait,
        n_variants = length(pair$ids), n = if (trait == "eqtl") n_eqtl else n_gwas,
        converged = isTRUE(f$converged), niter = f$niter,
        n_credible_sets = length(f$sets$cs), max_pip = max(f$pip),
        top_variant = names(f$pip)[which.max(f$pip)], final_elbo = tail(f$elbo, 1L),
        residual_variance = f$sigma2,
        finite_effective_rank = if (is.null(finite)) NA_real_ else finite$effective_rank,
        finite_rank_over_B = if (is.null(finite)) NA_real_ else finite$r_over_B,
        finite_max_penalty = if (is.null(finite)) NA_real_ else max(finite$per_variable_penalty),
        finite_sensitivity_flag = if (is.null(finite)) NA else finite$R_sensitivity_flag,
        finite_reliability_flag = if (is.null(finite)) NA else finite$R_reliability_flag,
        gwas_N_min = if (trait == "gwas") min(pair$gwas$N) else NA_real_,
        gwas_N_max = if (trait == "gwas") max(pair$gwas$N) else NA_real_
      )
    }
    fwrite(rbindlist(timing_rows), file.path(out_dir, "timings.tsv"), sep = "\t")
    message(format(Sys.time()), " | completed ", key)
  }
}

## write consolidated diagnostics only after every requested comparison finishes.
fwrite(rbindlist(fit_rows, fill = TRUE), file.path(out_dir, "fit_summary.tsv"), sep = "\t")
fwrite(rbindlist(fit_compare_rows, fill = TRUE), file.path(out_dir, "fit_route_comparison.tsv"), sep = "\t")
fwrite(rbindlist(coloc_rows, fill = TRUE), file.path(out_dir, "coloc_susie_summary.tsv"), sep = "\t")
fwrite(rbindlist(coloc_compare_rows, fill = TRUE), file.path(out_dir, "coloc_route_comparison.tsv"), sep = "\t")
fwrite(rbindlist(ld_rows, fill = TRUE), file.path(out_dir, "ld_diagnostics.tsv"), sep = "\t")
fwrite(rbindlist(exclusion_rows, fill = TRUE), file.path(out_dir, "exclusion_summary.tsv"), sep = "\t")
fwrite(rbindlist(timing_rows), file.path(out_dir, "timings.tsv"), sep = "\t")

## record approved input provenance and exact package installations.
input_paths <- c(
  paste0(plink_prefix, c(".pgen", ".pvar", ".psam")), ref_gt_file, ref_sample_file,
  vapply(disorders, gwas_bcf_path, character(1L), repo_root = repo_root)
)
q_paths <- file.path(tensorqtl_q_dir, sprintf("%s.tensorqtl_Q.tsv.gz", contexts))
if (nzchar(tensorqtl_q_dir)) input_paths <- c(input_paths, q_paths[file.exists(q_paths)])
provenance <- data.table(
  path = normalizePath(input_paths, mustWork = TRUE),
  bytes = file.info(input_paths)$size,
  sha256 = vapply(input_paths, sha256_file, character(1L))
)
fwrite(provenance, file.path(out_dir, "input_provenance.tsv"), sep = "\t")

packages <- c("coloc", "susieR", "data.table", "arrow", "pgenlibr")
package_provenance <- rbindlist(lapply(packages, function(pkg) {
  desc <- packageDescription(pkg)
  remote_sha <- if (is.null(desc$RemoteSha)) NA_character_ else desc$RemoteSha
  data.table(
    package = pkg, version = as.character(packageVersion(pkg)),
    library_path = find.package(pkg), remote_sha = remote_sha
  )
}), fill = TRUE)
fwrite(package_provenance, file.path(out_dir, "package_provenance.tsv"), sep = "\t")

parameters <- data.table(
  parameter = c("gene_id", "gene_name", "contexts", "disorders", "SI_min", "L",
                "scaled_prior_variance", "estimate_prior_variance",
                "estimate_residual_variance_eQTL", "estimate_residual_variance_GWAS",
                "R_finite_eQTL", "R_finite_GWAS", "coverage", "min_abs_corr",
                "max_iter", "z_method", "p1", "p2", "p12", "eqtl_LD",
                "eqtl_LD_sensitivities"),
  value = c(gene_id, gene_name, paste(contexts, collapse = ","), paste(disorders, collapse = ","),
            si_min, 10, 0.2, TRUE, TRUE, FALSE, FALSE, length(ref_samples),
            0.95, 0.5, 1000, "wald",
            1e-4, 1e-4, 1e-5, "rank-aware covariate-adjusted PGEN hardcalls",
            "raw donor LD; tensorQTL CPU-PyTorch QR LD")
)
fwrite(parameters, file.path(out_dir, "parameters.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Pilot outputs: ", out_dir)
