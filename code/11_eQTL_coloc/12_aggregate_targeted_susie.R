#!/usr/bin/env Rscript

## aggregate 95 checkpoints plus the reviewed MAPK3 pilot without rerunning models.

pilot_lib <- paste0(
    "/home/gpertea/work/R/spatialDLPFC_mdd_bpd/processed-data/11_eQTL_coloc/",
    "seurat/coloc/susie_exploratory/Rlib"
)
.libPaths(c(pilot_lib, .libPaths()))
suppressPackageStartupMessages({
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: 12_aggregate_targeted_susie.R BASE")
base <- normalizePath(args[[1L]], mustWork = TRUE)
aggregate_dir <- file.path(base, "aggregate")
dir.create(aggregate_dir, recursive = TRUE, showWarnings = FALSE)
manifest <- fread(file.path(base, "target_manifest.tsv"))
result_files <- file.path(base, "results", paste0(manifest$target_id, ".rds"))

mapk3_dir <- paste0(
    "/home/gpertea/work/R/spatialDLPFC_mdd_bpd/processed-data/11_eQTL_coloc/",
    "seurat/coloc/susie_targeted_no23andMe/MAPK3_pilot_20260819"
)
mapk3_fit_summary <- fread(file.path(mapk3_dir, "fit_summary.tsv"))
mapk3_ld <- fread(file.path(mapk3_dir, "ld_diagnostics.tsv"))
mapk3_exclusions <- fread(file.path(mapk3_dir, "exclusion_summary.tsv"))
mapk3_timings <- fread(file.path(mapk3_dir, "timings.tsv"))
mapk3_sensitivity <- fread(file.path(mapk3_dir, "coloc_p12_sensitivity.tsv"))

metadata_rows <- list()
fit_rows <- list()
ld_rows <- list()
exclusion_rows <- list()
coloc_rows <- list()
sensitivity_rows <- list()
comparison_rows <- list()
timing_rows <- list()
warning_rows <- list()
cs_rows <- list()
variant_audit_rows <- list()

## normalize one coloc summary and preserve no-credible-set outcomes explicitly
collect_coloc <- function(x, target_id, source) {
    s <- as.data.table(x$summary)
    if (!nrow(s)) {
        return(data.table(
            target_id = target_id, source = source,
            status = "no_credible_set_pair"
        ))
    }
    s[, `:=`(target_id = target_id, source = source, status = "credible_set_pair")]
    setcolorder(s, c("target_id", "source", "status"))
    s
}

## export credible-set membership and purity from compact SuSiE fits
collect_cs <- function(fit, target_id, trait) {
    sets <- fit$sets$cs
    if (is.null(sets) || !length(sets)) return(data.table())
    rbindlist(lapply(seq_along(sets), function(i) {
        idx <- sets[[i]]
        purity <- fit$sets$purity
        cs_label <- if (length(names(sets))) names(sets)[[i]] else paste0("L", i)
        data.table(
            target_id = target_id, trait = trait, cs = cs_label,
            cs_index = i, variant_id = names(fit$pip)[idx], pip = fit$pip[idx],
            cs_coverage = fit$sets$coverage[[i]],
            purity_min_abs_corr = if (is.null(purity)) NA_real_ else purity[i, "min.abs.corr"],
            purity_mean_abs_corr = if (is.null(purity)) NA_real_ else purity[i, "mean.abs.corr"],
            purity_median_abs_corr = if (is.null(purity)) NA_real_ else purity[i, "median.abs.corr"]
        )
    }), fill = TRUE)
}

for (i in seq_len(nrow(manifest))) {
    t <- manifest[i]
    id <- t$target_id
    if (!isTRUE(t$mapk3_pilot_complete)) {
        path <- result_files[[i]]
        if (!file.exists(path)) stop("missing checkpoint: ", id)
        x <- readRDS(path)
        metadata_rows[[id]] <- copy(x$metadata)
        fit_rows[[id]] <- copy(x$fit_summary)[, target_id := id]
        ld_rows[[id]] <- copy(x$ld_diagnostics)[, target_id := id]
        exclusion_rows[[id]] <- copy(x$exclusions)[, target_id := id]
        variant_audit_rows[[id]] <- copy(x$variant_audit)[, target_id := id]
        coloc_rows[[id]] <- collect_coloc(x$coloc_primary, id, "all_target_checkpoint")
        sensitivity_rows[[id]] <- copy(x$coloc_sensitivity)[, target_id := id]
        comparison_rows[[id]] <- copy(x$fit_comparison)[, target_id := id]
        timing_rows[[id]] <- copy(x$timings)[, target_id := id]
        if (length(x$warnings)) {
            warning_rows[[id]] <- data.table(target_id = id, warning = x$warnings)
        }
        cs_rows[[paste0(id, "_eqtl")]] <- collect_cs(x$fit_eqtl, id, "eqtl")
        cs_rows[[paste0(id, "_gwas")]] <- collect_cs(x$fit_gwas, id, "gwas")
    } else {
        key <- paste(t$disorder, t$dataset_id, sep = "__")
        x <- readRDS(file.path(mapk3_dir, "fits", paste0(key, ".rds")))
        metadata_rows[[id]] <- data.table(
            target_id = id, disorder = t$disorder, context = t$context,
            dataset_id = t$dataset_id, gene_id = t$gene_id, gene_name = t$gene_name,
            chr = t$chr, n_eqtl = unique(mapk3_fit_summary[comparison_key == key & trait == "eqtl", n]),
            covariate_rank = NA_integer_, n_reference = 503L,
            n_variants = length(x$ids), peak_rss_kb = NA_real_, completed = NA_character_
        )
        fit_rows[[id]] <- copy(mapk3_fit_summary[comparison_key == key])[, target_id := id]
        ld_rows[[id]] <- copy(mapk3_ld[comparison_key == key])[, target_id := id]
        exclusion_rows[[id]] <- copy(mapk3_exclusions[
            disorder == t$disorder & context == t$dataset_id
        ])[, target_id := id]
        variant_audit_rows[[id]] <- fread(file.path(
            mapk3_dir, "diagnostics", paste0(key, "_variant_audit.tsv.gz")
        ))[, target_id := id]
        coloc_rows[[id]] <- collect_coloc(x$coloc_prefit_primary, id, "MAPK3_reviewed_pilot")
        sensitivity_rows[[id]] <- copy(mapk3_sensitivity[
            disorder == t$disorder & context == t$dataset_id
        ])[, target_id := id]
        comparison_rows[[id]] <- data.table(
            target_id = id, comparison = "eqtl_rank_adjusted_vs_raw_LD",
            max_abs_pip_diff = max(abs(x$fit_eqtl_wrapper$pip - x$fit_eqtl_raw_ld$pip)),
            cs_identical = identical(x$fit_eqtl_wrapper$sets$cs, x$fit_eqtl_raw_ld$sets$cs),
            n_cs_primary = length(x$fit_eqtl_wrapper$sets$cs),
            n_cs_sensitivity = length(x$fit_eqtl_raw_ld$sets$cs)
        )
        timing_rows[[id]] <- copy(mapk3_timings[grepl(key, comparison_key, fixed = TRUE)])[
            , target_id := id
        ]
        cs_rows[[paste0(id, "_eqtl")]] <- collect_cs(x$fit_eqtl_wrapper, id, "eqtl")
        cs_rows[[paste0(id, "_gwas")]] <- collect_cs(x$fit_gwas_wrapper, id, "gwas")
    }
}

metadata <- rbindlist(metadata_rows, fill = TRUE)
fits <- rbindlist(fit_rows, fill = TRUE)
ld <- rbindlist(ld_rows, fill = TRUE)
exclusions <- rbindlist(exclusion_rows, fill = TRUE)
coloc <- rbindlist(coloc_rows, fill = TRUE)
sensitivity <- rbindlist(sensitivity_rows, fill = TRUE)
comparisons <- rbindlist(comparison_rows, fill = TRUE)
timings <- rbindlist(timing_rows, fill = TRUE)
warnings <- rbindlist(warning_rows, fill = TRUE)
credible_sets <- rbindlist(cs_rows, fill = TRUE)
variant_audit <- rbindlist(variant_audit_rows, fill = TRUE)

## summarize all credible-set pairs, retaining mixed H3/H4 signals within one target
pair_summary <- coloc[, .(
    n_coloc_pairs = sum(status == "credible_set_pair"),
    susie_max_H3 = if (any(status == "credible_set_pair")) max(PP.H3.abf, na.rm = TRUE) else NA_real_,
    susie_max_H4 = if (any(status == "credible_set_pair")) max(PP.H4.abf, na.rm = TRUE) else NA_real_,
    n_H3_ge_0.8 = sum(PP.H3.abf >= 0.8, na.rm = TRUE),
    n_H4_ge_0.8 = sum(PP.H4.abf >= 0.8, na.rm = TRUE)
), by = target_id]
fit_wide <- dcast(
    fits[, .(target_id, trait, converged, niter, n_credible_sets, max_pip, top_variant)],
    target_id ~ trait, value.var = c("converged", "niter", "n_credible_sets", "max_pip", "top_variant")
)

target_summary <- merge(manifest, metadata, by = intersect(names(manifest), names(metadata)), all.x = TRUE)
target_summary <- merge(target_summary, fit_wide, by = "target_id", all.x = TRUE)
target_summary <- merge(target_summary, pair_summary, by = "target_id", all.x = TRUE)
target_summary[, abf_H3_fraction := PP3 / PP34]
target_summary[, susie_interpretation := fcase(
    is.na(n_coloc_pairs) | n_coloc_pairs == 0L, "no_pair_no_resolvable_eqtl_gwas_cs",
    n_H4_ge_0.8 > 0L & n_H3_ge_0.8 > 0L, "mixed_shared_and_distinct_signal_pairs",
    n_H4_ge_0.8 > 0L, "shared_signal_pair_supported",
    n_H3_ge_0.8 > 0L, "distinct_signal_pair_supported",
    default = "credible_set_pairs_ambiguous"
)]
setorder(target_summary, target_id)

## parse external GNU time/status logs for full-process resource accounting
resource_rows <- lapply(manifest[mapk3_pilot_complete == FALSE, target_id], function(id) {
    err <- file.path(base, "logs", paste0(id, ".err"))
    status_file <- file.path(base, "logs", "status", paste0(id, ".tsv"))
    lines <- if (file.exists(err)) readLines(err, warn = FALSE) else character()
    status <- if (file.exists(status_file)) {
        fread(status_file, header = FALSE, sep = "\t", fill = TRUE)
    } else data.table()
    extract_num <- function(pattern) {
        x <- grep(pattern, lines, value = TRUE)
        if (!length(x)) return(NA_real_)
        value <- trimws(sub("^[^:]+:", "", tail(x, 1L)))
        as.numeric(sub("[^0-9.eE+-].*$", "", value))
    }
    data.table(
        target_id = id,
        exit_status = if (nrow(status)) status[[2L]][[1L]] else NA_integer_,
        user_seconds = extract_num("User time \\(seconds\\)"),
        system_seconds = extract_num("System time \\(seconds\\)"),
        cpu_percent = extract_num("Percent of CPU"),
        peak_rss_kb_external = extract_num("Maximum resident set size"),
        start = if (ncol(status) >= 3L) status[[3L]][[1L]] else NA_character_,
        end = if (ncol(status) >= 4L) status[[4L]][[1L]] else NA_character_
    )
})
resources <- rbindlist(resource_rows, fill = TRUE)

## write plain tables suitable for review without opening any fit objects
tables <- list(
    target_summary = target_summary, coloc_pairs = coloc,
    fit_summary = fits, ld_diagnostics = ld, exclusion_summary = exclusions,
    fit_sensitivity = comparisons, p12_sensitivity = sensitivity,
    credible_set_members = credible_sets, timings = timings,
    warnings = warnings, resources = resources, variant_audit = variant_audit
)
for (name in names(tables)) {
    fwrite(tables[[name]], file.path(aggregate_dir, paste0(name, ".tsv.gz")), sep = "\t")
}

## exact cardinality and convergence checks are terminal audit requirements
audit <- data.table(
    check = c(
        "manifest_rows", "target_summary_rows", "unique_target_ids", "non_mapk3_checkpoints",
        "trait_fit_rows", "nonconverged_fits", "failed_processes", "missing_exclusion_ledgers"
    ),
    value = c(
        nrow(manifest), nrow(target_summary), uniqueN(target_summary$target_id),
        sum(file.exists(result_files[!manifest$mapk3_pilot_complete])),
        nrow(fits), sum(!fits$converged), sum(resources$exit_status != 0L, na.rm = TRUE),
        sum(!manifest$target_id %chin% exclusions$target_id)
    )
)
fwrite(audit, file.path(aggregate_dir, "completion_audit.tsv"), sep = "\t")
