#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(SeqArray)
})

## extract exact target variants and approved samples from a fresh all-chromosome GDS
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: 09_extract_1000g_target_reference.R BASE")
base <- normalizePath(args[[1L]], mustWork = TRUE)
ref_dir <- file.path(base, "inputs", "reference")
audit_dir <- file.path(base, "audit", "reference")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

gds_url <- paste0(
    "https://gds-stat.s3.amazonaws.com/download/1000g/2022/",
    "1kGP_high_coverage_Illumina.allchr.filtered.SNV_INDEL_SV_phased_panel.gds"
)
gds_file <- file.path(ref_dir, basename(gds_url))
sample_file <- paste0(
    "/home/gpertea/work/R/spatialDLPFC_mdd_bpd/processed-data/11_eQTL_coloc/",
    "seurat/coloc/susie_exploratory/reference/1000G_phase3_EUR_2504-unrelated.samples"
)
target_file <- file.path(base, "inputs", "target_pvar.tsv.gz")
stopifnot(file.exists(gds_file), file.exists(sample_file), file.exists(target_file))

samples <- scan(sample_file, what = "", quiet = TRUE)
target <- fread(target_file)
stopifnot(length(samples) == 503L, !anyDuplicated(samples), !anyDuplicated(target$variant_id))
target[, chr_code := sub("^chr", "", CHROM)]
target[, chr_num := as.integer(chr_code)]

gds <- seqOpen(gds_file)
on.exit(seqClose(gds), add = TRUE)
all_samples <- seqGetData(gds, "sample.id")
stopifnot(all(samples %chin% all_samples), length(all_samples) == 3202L)

## preserve the approved sample order and iterate by chromosome to bound memory
seqSetFilter(gds, NULL, sample.id = samples, verbose = FALSE)
audit <- vector("list", length(unique(target$chr_num)))
chromosomes <- unique(target[order(chr_num), CHROM])

for (i in seq_along(chromosomes)) {
    chr <- chromosomes[[i]]
    x <- target[CHROM == chr]
    seqResetFilter(gds, sample = FALSE, variant = TRUE, verbose = FALSE)
    idx <- seqSetFilterPos(
        gds, x$chr_code, x$POS, x$REF, x$ALT,
        intersect = TRUE, ret.idx = TRUE, verbose = FALSE
    )

    ## GDS $dosage counts reference alleles; convert to ALT allele count
    pos <- seqGetData(gds, "position")
    allele <- tstrsplit(seqGetData(gds, "allele"), ",", fixed = TRUE)
    dosage <- 2 - seqGetData(gds, "$dosage")
    ids <- paste0(chr, ":", pos, ":", allele[[1L]], ":", allele[[2L]])
    stopifnot(ncol(dosage) == length(ids), all(ids %chin% x$variant_id))

    ## write variants by genomic order with numeric ALT dosages for fast LD loading
    out <- as.data.table(t(dosage))
    setnames(out, samples)
    out[, variant_id := ids]
    setcolorder(out, c("variant_id", samples))
    out_file <- file.path(ref_dir, paste0(chr, ".1000G_highcov.EUR_unrelated.GT.tsv.gz"))
    fwrite(out, out_file, sep = "\t", na = "NA", compress = "gzip")

    audit[[i]] <- data.table(
        chr = chr,
        n_requested = nrow(x),
        n_exact_match = length(ids),
        n_absent = sum(is.na(idx)),
        n_samples = nrow(dosage),
        n_missing_genotypes = sum(is.na(dosage)),
        gt_file = out_file,
        gt_md5 = unname(tools::md5sum(out_file))
    )
    message(chr, ": ", length(ids), "/", nrow(x), " exact variants")
}

audit <- rbindlist(audit)
fwrite(audit, file.path(audit_dir, "reference_extraction.tsv"), sep = "\t")

## record immutable source and package provenance separately from derived products
prov <- data.table(
    field = c(
        "source_url", "source_description", "source_size_bytes", "source_md5",
        "gds_samples", "gds_variants", "approved_samples", "SeqArray", "gdsfmt", "R"
    ),
    value = c(
        gds_url,
        "UW SeqArray conversion of official 20220422 1kGP 30x phased panel",
        file.info(gds_file)$size,
        unname(tools::md5sum(gds_file)),
        length(all_samples), 73554796L, length(samples),
        as.character(packageVersion("SeqArray")),
        as.character(packageVersion("gdsfmt")),
        as.character(getRversion())
    )
)
fwrite(prov, file.path(audit_dir, "reference_provenance.tsv"), sep = "\t")
