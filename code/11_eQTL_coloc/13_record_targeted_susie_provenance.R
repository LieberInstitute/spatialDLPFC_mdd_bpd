#!/usr/bin/env Rscript

## record exact files, package builds, host, and authoritative source URLs.

pilot_lib <- paste0(
    "/home/gpertea/work/R/spatialDLPFC_mdd_bpd/processed-data/11_eQTL_coloc/",
    "seurat/coloc/susie_exploratory/Rlib"
)
.libPaths(c(pilot_lib, .libPaths()))
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: 13_record_targeted_susie_provenance.R BASE")
base <- normalizePath(args[[1L]], mustWork = TRUE)
audit_dir <- file.path(base, "audit")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
manifest <- fread(file.path(base, "target_manifest.tsv"))
gwas_manifest <- fread(file.path(base, "gwas_manifest.tsv"))
ref_manifest <- fread(file.path(base, "reference_request_manifest.tsv"))

repo_root <- "/home/gpertea/work/R/spatialDLPFC_mdd_bpd"
processed <- file.path(repo_root, "processed-data", "11_eQTL_coloc", "seurat")
plink_prefix <- file.path(repo_root, "processed-data", "00_genotypes", "plink2", "merged_maf05")
sample_file <- file.path(
    processed, "coloc", "susie_exploratory", "reference",
    "1000G_phase3_EUR_2504-unrelated.samples"
)

## include every direct model input while excluding checkpoints and old caches
paths <- rbindlist(list(
    data.table(role = "approved_GWAS_BCF", path = gwas_manifest$approved_bcf),
    data.table(role = "fresh_GWAS_approved_extract", path = gwas_manifest$approved_file),
    data.table(role = "fresh_GWAS_harmonization", path = gwas_manifest$harmonization_file),
    data.table(role = "target_nominal_tensorQTL", path = unique(manifest$nominal_file)),
    data.table(role = "target_covariates", path = file.path(
        processed, "tqtl_in", paste0(unique(manifest$dataset_id), ".gene.covars.txt")
    )),
    data.table(role = "PGEN_hardcalls", path = paste0(plink_prefix, c(".pgen", ".pvar", ".psam"))),
    data.table(role = "fresh_1000G_GDS", path = file.path(
        base, "inputs", "reference",
        "1kGP_high_coverage_Illumina.allchr.filtered.SNV_INDEL_SV_phased_panel.gds"
    )),
    data.table(role = "fresh_1000G_target_GT", path = ref_manifest$gt_file),
    data.table(role = "approved_1000G_samples", path = sample_file),
    data.table(role = "target_PVAR", path = file.path(base, "inputs", "target_pvar.tsv.gz")),
    data.table(role = "target_manifest", path = file.path(base, "target_manifest.tsv")),
    data.table(role = "final_coloc_pass_source", path = file.path(processed, "final", "coloc_results.xlsx"))
), use.names = TRUE)
paths <- unique(paths[file.exists(path)])
paths[, path := normalizePath(path, mustWork = TRUE)]

sha256_file <- function(path) {
    x <- system2("sha256sum", shQuote(path), stdout = TRUE)
    sub("[[:space:]].*$", "", x[[1L]])
}
message("hashing ", nrow(paths), " direct input files")
paths[, `:=`(
    bytes = file.info(path)$size,
    modified = format(file.info(path)$mtime, "%Y-%m-%d %H:%M:%S %z"),
    sha256 = vapply(path, sha256_file, character(1L))
)]
setorder(paths, role, path)
fwrite(paths, file.path(audit_dir, "input_provenance.tsv.gz"), sep = "\t")

packages <- c("coloc", "susieR", "data.table", "pgenlibr", "SeqArray", "gdsfmt")
package_provenance <- rbindlist(lapply(packages, function(pkg) {
    desc <- packageDescription(pkg)
    data.table(
        package = pkg, version = as.character(packageVersion(pkg)),
        library_path = find.package(pkg),
        remote_sha = if (is.null(desc$RemoteSha)) NA_character_ else desc$RemoteSha
    )
}), fill = TRUE)
fwrite(package_provenance, file.path(audit_dir, "package_provenance.tsv"), sep = "\t")

sources <- data.table(
    resource = c(
        "official_1000G_30x_collection", "UW_SeqArray_catalog", "fresh_GDS_object",
        "susieR_RSS", "susieR_RSS_diagnostics", "susieR_LD_mismatch",
        "coloc_SuSiE", "tensorQTL_source", "TOP-LD_publication", "TOP-LD_API_client"
    ),
    url = c(
        "https://internationalgenome.org/data-portal/data-collections/1000genomes_30x/",
        "https://gds-stat.s3.amazonaws.com/download/1000g/index.html",
        paste0("https://gds-stat.s3.amazonaws.com/download/1000g/2022/",
               "1kGP_high_coverage_Illumina.allchr.filtered.SNV_INDEL_SV_phased_panel.gds"),
        "https://stephenslab.github.io/susieR/reference/susie_rss.html",
        "https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html",
        "https://stephenslab.github.io/susieR/articles/rss_mismatch.html",
        "https://chr1swallace.github.io/coloc/articles/a06_SuSiE.html",
        "https://github.com/broadinstitute/tensorqtl/blob/master/tensorqtl/core.py",
        "https://doi.org/10.1016/j.ajhg.2022.08.006",
        "https://github.com/linnabrown/topld_api"
    )
)
fwrite(sources, file.path(audit_dir, "source_documentation.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(audit_dir, "sessionInfo.txt"))
writeLines(capture.output(Sys.info()), file.path(audit_dir, "host_info.txt"))
