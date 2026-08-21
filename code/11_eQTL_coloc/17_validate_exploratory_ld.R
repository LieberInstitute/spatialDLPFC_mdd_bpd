#!/usr/bin/env Rscript

## report signed donor and reference LD for predeclared variants and SuSiE sets.

suppressPackageStartupMessages({
  library(data.table)
  library(pgenlibr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: 17_validate_exploratory_ld.R BASE")
base <- normalizePath(args[[1L]], mustWork = TRUE)
repo_root <- "/home/gpertea/work/R/spatialDLPFC_mdd_bpd"
processed <- file.path(repo_root, "processed-data", "11_eQTL_coloc", "seurat")
out_dir <- file.path(base, "ld_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
manifest <- fread(file.path(base, "target_manifest.tsv"))
pvar <- fread(file.path(base, "inputs", "target_pvar.tsv.gz"))
sig <- fread(file.path(processed, "tables", "map_significant_pairs.csv.gz"))
sig <- sig[split == "all"]

plink_prefix <- file.path(repo_root, "processed-data", "00_genotypes", "plink2", "merged_maf05")
psam <- fread(paste0(plink_prefix, ".psam"))
sample_col <- intersect(c("#IID", "IID"), names(psam))[[1L]]
pvar_object <- NewPvar(paste0(plink_prefix, ".pvar"))
pgen <- NewPgen(paste0(plink_prefix, ".pgen"), pvar = pvar_object)
on.exit(try(ClosePgen(pgen), silent = TRUE), add = TRUE)

pairwise_rows <- list()
variant_rows <- list()
summary_rows <- list()

fit_cs <- function(fit) {
  cs <- fit$sets$cs
  if (is.null(cs) || !length(cs)) return(character())
  unique(unlist(lapply(cs, function(i) names(fit$pip)[i]), use.names = FALSE))
}

plot_ld <- function(R, target_id, source) {
  n <- ncol(R)
  labels <- if (n <= 50L) sub("^chr[^:]+:", "", colnames(R)) else rep("", n)
  pdf(file.path(out_dir, sprintf("%s_%s_signed_r.pdf", target_id, source)),
      width = if (n <= 50L) 10 else 8, height = if (n <= 50L) 10 else 8)
  on.exit(dev.off())
  cols <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(201)
  par(mar = if (n <= 50L) c(9, 9, 3, 2) else c(3, 3, 3, 2))
  image(seq_len(n), seq_len(n), t(R[n:1, , drop = FALSE]), zlim = c(-1, 1),
        col = cols, axes = FALSE, xlab = "", ylab = "",
        main = paste(target_id, source, "signed r"))
  if (n <= 50L) {
    axis(1, at = seq_len(n), labels = labels, las = 2, cex.axis = 0.45)
    axis(2, at = seq_len(n), labels = rev(labels), las = 2, cex.axis = 0.45)
  }
  box()
}

for (i in seq_len(nrow(manifest))) {
  target <- manifest[i]
  none <- readRDS(file.path(base, "results", paste0(target$target_id, ".rds")))
  eb <- readRDS(file.path(base, "results_eb", paste0(target$target_id, ".rds")))
  included <- none$variant_audit[is.na(exclusion_reason), variant_id]
  reported <- sig[dataset_id == target$dataset_id & gene_id == target$gene_id, variant_id]
  roles <- list(
    reported_eqtl = reported,
    eqtl_cs = fit_cs(none$fit_eqtl),
    gwas_cs_none = fit_cs(none$fit_gwas),
    gwas_cs_eb = fit_cs(eb$fit_gwas),
    eqtl_top = names(none$fit_eqtl$pip)[which.max(none$fit_eqtl$pip)],
    gwas_top_none = names(none$fit_gwas$pip)[which.max(none$fit_gwas$pip)],
    gwas_top_eb = names(eb$fit_gwas$pip)[which.max(eb$fit_gwas$pip)]
  )
  focus <- unique(unlist(roles, use.names = FALSE))
  focus <- focus[focus %chin% included]
  if (length(focus) < 2L) next
  pv <- pvar[match(focus, variant_id)]
  setorder(pv, POS, var_idx)
  focus <- pv$variant_id

  ## calculate raw donor LD and the covariate-adjusted LD used by eQTL SuSiE.
  X_all <- ReadIntList(pgen, pv$var_idx)
  rownames(X_all) <- as.character(psam[[sample_col]])
  colnames(X_all) <- focus
  cov_raw <- fread(file.path(processed, "tqtl_in", sprintf("%s.gene.covars.txt", target$dataset_id)))
  C <- t(as.matrix(cov_raw[, -1L]))
  storage.mode(C) <- "double"
  rownames(C) <- names(cov_raw)[-1L]
  X <- X_all[rownames(C), focus, drop = FALSE]
  for (j in which(colSums(is.na(X)) > 0L)) X[is.na(X[, j]), j] <- mean(X[, j], na.rm = TRUE)
  R_raw <- cor(X)
  R_adjusted <- cor(qr.resid(qr(C, tol = 1e-7, LAPACK = FALSE), X))

  ## calculate external-reference LD with the exact 503 approved EUR samples.
  ref <- fread(file.path(base, "inputs", "reference",
                         sprintf("%s.1000G_highcov.EUR_unrelated.GT.tsv.gz", target$chr)))
  ref_samples <- setdiff(names(ref), "variant_id")
  X_ref <- t(as.matrix(ref[match(focus, variant_id), ..ref_samples]))
  storage.mode(X_ref) <- "double"
  colnames(X_ref) <- focus
  for (j in which(colSums(is.na(X_ref)) > 0L)) {
    X_ref[is.na(X_ref[, j]), j] <- mean(X_ref[, j], na.rm = TRUE)
  }
  R_ref <- cor(X_ref)

  role_table <- data.table(variant_id = focus, pos = pv$POS)
  for (role in names(roles)) role_table[, (role) := variant_id %chin% roles[[role]]]
  role_table[, `:=`(target_id = target$target_id, gene_name = target$gene_name,
                    context = target$context)]
  variant_rows[[target$target_id]] <- role_table

  for (source in c("donor_raw", "donor_adjusted", "reference_503_EUR")) {
    R <- switch(source, donor_raw = R_raw, donor_adjusted = R_adjusted,
                reference_503_EUR = R_ref)
    idx <- which(upper.tri(R), arr.ind = TRUE)
    pairs <- data.table(
      target_id = target$target_id, gene_name = target$gene_name,
      context = target$context, source = source,
      variant1 = colnames(R)[idx[, 1L]], variant2 = colnames(R)[idx[, 2L]],
      r = R[idx], r2 = R[idx]^2,
      distance_bp = abs(pv$POS[idx[, 1L]] - pv$POS[idx[, 2L]])
    )
    pairwise_rows[[paste(target$target_id, source)]] <- pairs
    summary_rows[[paste(target$target_id, source)]] <- pairs[, .(
      n_focus_variants = ncol(R), n_pairs = .N,
      n_abs_r_ge_0_8 = sum(abs(r) >= 0.8), n_abs_r_ge_0_99 = sum(abs(r) >= 0.99),
      n_exact_abs_r_1 = sum(abs(abs(r) - 1) < 1e-12),
      min_r = min(r), max_r = max(r), median_abs_r = median(abs(r))
    ), by = .(target_id, gene_name, context, source)]
    plot_ld(R, target$target_id, source)
  }
}

fwrite(rbindlist(variant_rows, fill = TRUE), file.path(out_dir, "focus_variants.tsv"), sep = "\t")
fwrite(rbindlist(pairwise_rows, fill = TRUE), file.path(out_dir, "focus_pairwise_ld.tsv.gz"), sep = "\t")
fwrite(rbindlist(summary_rows, fill = TRUE), file.path(out_dir, "focus_ld_summary.tsv"), sep = "\t")
message("completed exploratory signed-LD validation")
