#!/usr/bin/env Rscript

## aggregate isolated baseline and EB checkpoints into small review tables.

pilot_lib <- paste0(
  "/home/gpertea/work/R/spatialDLPFC_mdd_bpd/processed-data/11_eQTL_coloc/",
  "seurat/coloc/susie_exploratory/Rlib"
)
.libPaths(c(pilot_lib, .libPaths()))
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: 16_aggregate_susie_next.R BASE")
base <- normalizePath(args[[1L]], mustWork = TRUE)
out_dir <- file.path(base, "summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
manifest <- fread(file.path(base, "target_manifest.tsv"))

canonical_cs <- function(fit, trait, mode, target) {
  cs <- fit$sets$cs
  if (is.null(cs) || !length(cs)) return(data.table())
  rbindlist(lapply(seq_along(cs), function(i) {
    idx <- cs[[i]]
    data.table(
      target_id = target$target_id, disorder = target$disorder,
      context = target$context, gene_id = target$gene_id,
      gene_name = target$gene_name, trait = trait, mode = mode,
      cs = names(cs)[i], variant_id = names(fit$pip)[idx], pip = fit$pip[idx],
      cs_size = length(idx), cs_min_abs_corr = fit$sets$purity$min.abs.corr[i],
      cs_mean_abs_corr = fit$sets$purity$mean.abs.corr[i],
      cs_median_abs_corr = fit$sets$purity$median.abs.corr[i]
    )
  }))
}

collect_h4 <- function(coloc_result, target_id, mode) {
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
      target_id = target_id, mode = mode, signal_pair = i,
      hit1 = summary$hit1[[i]], hit2 = summary$hit2[[i]],
      pair_H3 = summary$PP.H3.abf[[i]], pair_H4 = summary$PP.H4.abf[[i]]
    )]
    x
  }), fill = TRUE)
}

rows <- list(fits = list(), comparisons = list(), coloc = list(),
             cs = list(), h4_variants = list(), exclusions = list(), ld = list(), timings = list(),
             warnings = list(), metadata = list())

for (mode in c("none", "eb")) {
  result_dir <- file.path(base, if (mode == "none") "results" else "results_eb")
  files <- file.path(result_dir, paste0(manifest$target_id, ".rds"))
  if (!all(file.exists(files))) stop("Incomplete ", mode, " checkpoint set")
  for (i in seq_len(nrow(manifest))) {
    target <- manifest[i]
    x <- readRDS(files[[i]])
    key <- paste(target$target_id, mode, sep = "__")
    fs <- copy(x$fit_summary)
    fs[, `:=`(target_id = target$target_id, mode = mode)]
    rows$fits[[key]] <- fs
    dc <- copy(x$direct_wrapper_comparison)
    dc[, `:=`(target_id = target$target_id, mode = mode)]
    rows$comparisons[[key]] <- dc
    co <- as.data.table(x$coloc_primary$summary)
    if (!nrow(co)) co <- data.table(status = "no_credible_set_pair")
    co[, `:=`(target_id = target$target_id, mode = mode)]
    rows$coloc[[key]] <- co
    rows$h4_variants[[key]] <- collect_h4(x$coloc_primary, target$target_id, mode)
    rows$cs[[paste0(key, "__eqtl")]] <- canonical_cs(x$fit_eqtl, "eqtl", mode, target)
    rows$cs[[paste0(key, "__gwas")]] <- canonical_cs(x$fit_gwas, "gwas", mode, target)
    ex <- copy(x$exclusions)
    ex[, `:=`(target_id = target$target_id, mode = mode)]
    rows$exclusions[[key]] <- ex
    ld <- copy(x$ld_diagnostics)
    ld[, `:=`(target_id = target$target_id, mode = mode)]
    rows$ld[[key]] <- ld
    tm <- copy(x$timings)
    tm[, `:=`(target_id = target$target_id, mode = mode)]
    rows$timings[[key]] <- tm
    rows$warnings[[key]] <- data.table(
      target_id = target$target_id, mode = mode, warning = x$warnings
    )
    md <- copy(x$metadata)
    md[, mode := mode]
    rows$metadata[[key]] <- md
  }
}

for (name in names(rows)) {
  value <- rbindlist(rows[[name]], fill = TRUE)
  fwrite(value, file.path(out_dir, paste0(name, ".tsv")), sep = "\t")
}

## compare baseline and EB at matched signal-pair labels where possible.
coloc <- rbindlist(rows$coloc, fill = TRUE)
id_cols <- intersect(c("target_id", "hit1", "hit2"), names(coloc))
metric_cols <- intersect(c("PP.H3.abf", "PP.H4.abf"), names(coloc))
if (length(metric_cols)) {
  wide <- dcast(coloc, as.formula(paste(paste(id_cols, collapse = " + "), "~ mode")),
                value.var = metric_cols)
  fwrite(wide, file.path(out_dir, "coloc_none_vs_eb.tsv"), sep = "\t")
}

message("aggregated ", nrow(manifest), " targets across baseline and EB")
