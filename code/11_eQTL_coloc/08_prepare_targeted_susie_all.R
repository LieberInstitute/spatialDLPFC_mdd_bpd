#!/usr/bin/env Rscript

## prepare approved targets without reading production coloc caches.
## optional arguments allow a reviewed exploratory manifest and isolated output.

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(dplyr)
  library(openxlsx)
})

repo_root <- normalizePath(getwd(), mustWork = TRUE)
if (!dir.exists(file.path(repo_root, ".git"))) stop("Run from repository root")
source(file.path(repo_root, "code", "11_eQTL_coloc", "utils.R"))
data.table::setDTthreads(2L)

base_dir <- file.path(repo_root, "processed-data", "11_eQTL_coloc", "seurat")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) {
  stop("usage: 08_prepare_targeted_susie_all.R [TARGET_MANIFEST.tsv OUT_DIR]")
}
if (length(args) == 1L) stop("TARGET_MANIFEST and OUT_DIR must be supplied together")
exploratory <- length(args) == 2L
out_dir <- if (exploratory) args[[2L]] else file.path(
  base_dir, "coloc", "susie_targeted_no23andMe", "all_targets_20260819"
)
nominal_dir <- file.path(out_dir, "inputs", "nominal")
gwas_dir <- file.path(out_dir, "inputs", "gwas")
ref_dir <- file.path(out_dir, "inputs", "reference")
dir.create(nominal_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gwas_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)

if (exploratory) {
  manifest_file <- normalizePath(args[[1L]], mustWork = TRUE)
  targets <- fread(manifest_file)
  required <- c("target_id", "disorder", "context", "dataset_id", "gene_id",
                "gene_name", "chr", "region_start", "region_end")
  if (!all(required %in% names(targets))) stop("Exploratory manifest is incomplete")
  if (anyDuplicated(targets$target_id)) stop("Exploratory target IDs are not unique")
  setnames(targets, c("region_start", "region_end"),
           c("declared_region_start", "declared_region_end"))
  targets[, mapk3_pilot_complete := FALSE]
} else {
  final_xlsx <- file.path(base_dir, "final", "coloc_results.xlsx")
  targets <- as.data.table(read.xlsx(final_xlsx, sheet = "coloc_pass"))
  if (nrow(targets) != 100L) stop("Expected exactly 100 coloc_pass rows")
  targets[, target_id := sprintf("target_%03d", .I)]
  targets[, dataset_id := unname(SEURAT_CONTEXT_TO_DATASET_ID[context])]
  if (anyNA(targets$dataset_id)) stop("Unmapped final-table context")
  targets[, chr := sub(":.*$", "", lead_snp)]
  targets[, mapk3_pilot_complete := gene_id == "ENSG00000102882"]
}
targets[, chr_num := as.integer(sub("^chr", "", chr))]
if (anyNA(targets$chr_num)) stop("Cannot parse target chromosome")

## scan each required parquet once and export only target genes.
gene_context <- unique(targets[, .(
  dataset_id, context, gene_id, gene_name, chr, chr_num
)])
gene_context[, nominal_file := file.path(
  nominal_dir, sprintf("%s__%s.tsv.gz", dataset_id, gene_id)
)]

nominal_cols <- c(
  "phenotype_id", "variant_id", "start_distance", "af", "ma_samples",
  "ma_count", "pval_nominal", "slope", "slope_se"
)
for (i in seq_len(nrow(unique(gene_context[, .(dataset_id, chr_num)])))) {
  group <- unique(gene_context[, .(dataset_id, chr_num)])[i]
  genes <- gene_context[dataset_id == group$dataset_id & chr_num == group$chr_num, gene_id]
  parquet <- file.path(
    base_dir, "tqtl_out",
    sprintf("%s.gene.cis_qtl_pairs.chr%d.parquet", group$dataset_id, group$chr_num)
  )
  if (!file.exists(parquet)) stop("Missing nominal parquet: ", parquet)
  message("reading ", basename(parquet), " for ", length(genes), " target genes")
  tab <- open_dataset(parquet) |>
    filter(phenotype_id %in% genes) |>
    select(all_of(nominal_cols)) |>
    collect()
  tab <- as.data.table(tab)
  for (gene in genes) {
    x <- tab[phenotype_id == gene]
    if (!nrow(x)) stop("Missing nominal target ", gene, " in ", parquet)
    setorder(x, variant_id)
    fwrite(x, gene_context[dataset_id == group$dataset_id & gene_id == gene, nominal_file], sep = "\t")
  }
}

## collect the exact variant universe and map it once to PGEN indexes/alleles.
variant_lists <- lapply(gene_context$nominal_file, function(path) fread(path, select = "variant_id")$variant_id)
names(variant_lists) <- paste(gene_context$dataset_id, gene_context$gene_id, sep = "__")
all_ids <- unique(unlist(variant_lists, use.names = FALSE))
plink_prefix <- file.path(repo_root, "processed-data", "00_genotypes", "plink2", "merged_maf05")
pvar <- read_plink2_variant_table(plink_prefix)
pvar_target <- pvar[data.table(variant_id = all_ids), on = "variant_id", nomatch = 0L]
if (nrow(pvar_target) != length(all_ids)) stop("Target nominal variants missing from PVAR")
setorder(pvar_target, var_idx)
pvar_file <- file.path(out_dir, "inputs", "target_pvar.tsv.gz")
fwrite(pvar_target, pvar_file, sep = "\t")

## define exact gene intervals and collapse overlaps only for input reuse.
interval_rows <- lapply(seq_len(nrow(gene_context)), function(i) {
  g <- gene_context[i]
  ids <- variant_lists[[paste(g$dataset_id, g$gene_id, sep = "__")]]
  pos <- pvar_target[match(ids, variant_id), POS]
  data.table(
    dataset_id = g$dataset_id, context = g$context, gene_id = g$gene_id,
    gene_name = g$gene_name, chr = g$chr, chr_num = g$chr_num,
    nominal_file = g$nominal_file, n_nominal = length(ids),
    region_start = min(pos), region_end = max(pos)
  )
})
gene_context <- rbindlist(interval_rows)

gene_intervals <- gene_context[, .(
  gene_name = gene_name[[1L]], chr = chr[[1L]], chr_num = chr_num[[1L]],
  region_start = min(region_start), region_end = max(region_end)
), by = gene_id]
setorder(gene_intervals, chr_num, region_start, region_end)
gene_intervals[, region_id := NA_character_]
region_rows <- list()
for (chrom in unique(gene_intervals$chr_num)) {
  idx <- which(gene_intervals$chr_num == chrom)
  region_no <- 0L
  current_end <- -Inf
  current_id <- ""
  for (j in idx) {
    if (gene_intervals$region_start[[j]] > current_end) {
      region_no <- region_no + 1L
      current_id <- sprintf("chr%02d_region%02d", chrom, region_no)
      current_end <- gene_intervals$region_end[[j]]
      region_rows[[current_id]] <- data.table(
        region_id = current_id, chr = gene_intervals$chr[[j]], chr_num = chrom,
        region_start = gene_intervals$region_start[[j]], region_end = current_end
      )
    } else {
      current_end <- max(current_end, gene_intervals$region_end[[j]])
      region_rows[[current_id]][, region_end := current_end]
    }
    gene_intervals$region_id[[j]] <- current_id
  }
}
regions <- rbindlist(region_rows)
gene_context <- merge(
  gene_context, gene_intervals[, .(gene_id, region_id)],
  by = "gene_id", all.x = TRUE
)

## attach reusable region and file definitions to every final comparison row.
targets <- merge(
  targets,
  gene_context[, .(dataset_id, gene_id, nominal_file, n_nominal,
                    region_start, region_end, region_id)],
  by = c("dataset_id", "gene_id"), all.x = TRUE
)
setorder(targets, target_id)
if (anyNA(targets$nominal_file)) stop("Incomplete target-to-nominal mapping")
if (exploratory && any(
  targets$declared_region_start != targets$region_start |
  targets$declared_region_end != targets$region_end
)) stop("Reviewed exploratory bounds differ from dense nominal bounds")

## query each approved GWAS BCF once at all target positions, without p filtering.
query_gwas <- function(disorder) {
  bcf <- gwas_bcf_path(disorder, repo_root = repo_root)
  needed_genes <- unique(targets[get("disorder") == disorder, gene_id])
  needed_ids <- unique(unlist(variant_lists[
    vapply(strsplit(names(variant_lists), "__", fixed = TRUE), `[[`, character(1L), 2L) %in% needed_genes
  ], use.names = FALSE))
  pv <- pvar_target[variant_id %in% needed_ids]
  bed <- file.path(gwas_dir, sprintf("%s_target_positions.bed", disorder))
  raw_file <- file.path(gwas_dir, sprintf("%s_query_raw.tsv.gz", disorder))
  all_file <- file.path(gwas_dir, sprintf("%s_harmonized_all.tsv.gz", disorder))
  keep_file <- file.path(gwas_dir, sprintf("%s_approved.tsv.gz", disorder))
  fwrite(unique(pv[, .(CHROM, start0 = POS - 1L, end1 = POS)]), bed,
         sep = "\t", col.names = FALSE)

  tmp <- tempfile(fileext = ".tsv")
  err <- tempfile(fileext = ".log")
  on.exit(unlink(c(tmp, err)), add = TRUE)
  status <- system2(
    "bcftools",
    c("query", "-R", shQuote(bed),
      "-f", shQuote("%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%ES\t%SE\t%LP\t%NE\t%NS\t%NC\t%SI]\n"),
      shQuote(bcf)), stdout = tmp, stderr = err
  )
  if (status != 0L) stop(paste(readLines(err, warn = FALSE), collapse = "\n"))
  raw <- fread(tmp, header = FALSE)
  setnames(raw, c("chr", "pos", "rsid", "a0", "a1", "beta", "beta_se",
                  "lp", "N", "ns", "ncas", "impinfo"))
  raw[, `:=`(
    chr = as.character(chr), pos = as.integer(pos), a0 = as.character(a0), a1 = as.character(a1),
    beta = suppressWarnings(as.numeric(beta)), beta_se = suppressWarnings(as.numeric(beta_se)),
    lp = suppressWarnings(as.numeric(lp)), N = suppressWarnings(as.numeric(N)),
    ncas = suppressWarnings(as.numeric(ncas)), impinfo = suppressWarnings(as.numeric(impinfo)),
    p = 10^(-suppressWarnings(as.numeric(lp)))
  )]
  fwrite(raw, raw_file, sep = "\t")

  exact <- merge(raw, pv, by.x = c("chr", "pos", "a0", "a1"),
                 by.y = c("CHROM", "POS", "REF", "ALT"), allow.cartesian = TRUE)
  exact[, `:=`(match_mode = "exact", match_rank = 1L)]
  swapped <- merge(raw, pv, by.x = c("chr", "pos", "a0", "a1"),
                   by.y = c("CHROM", "POS", "ALT", "REF"), allow.cartesian = TRUE)
  swapped[, `:=`(beta = -beta, match_mode = "swapped", match_rank = 2L)]
  matched <- rbindlist(list(exact, swapped), use.names = TRUE, fill = TRUE)
  setorder(matched, variant_id, match_rank, p)
  matched <- matched[!duplicated(variant_id)]
  matched[, exclusion_reason := fcase(
    !is.finite(impinfo), "missing_SI",
    impinfo < 0.8, "SI_below_0.8",
    !is.finite(beta) | !is.finite(beta_se) | beta_se <= 0, "invalid_beta_or_se",
    !is.finite(N) | N <= 0, "invalid_N",
    !is.finite(p) | p <= 0 | p > 1, "invalid_p",
    default = NA_character_
  )]
  fwrite(matched, all_file, sep = "\t")
  fwrite(matched[is.na(exclusion_reason)], keep_file, sep = "\t")
  data.table(
    disorder = disorder, approved_bcf = normalizePath(bcf),
    n_requested = length(needed_ids), n_raw = nrow(raw),
    n_allele_matched = nrow(matched), n_approved = sum(is.na(matched$exclusion_reason)),
    approved_file = keep_file, harmonization_file = all_file
  )
}

gwas_rows <- list()
for (disorder_value in unique(targets$disorder)) {
  message("extracting approved GWAS ", disorder_value)
  gwas_rows[[disorder_value]] <- query_gwas(disorder_value)
}
gwas_manifest <- rbindlist(gwas_rows)

## prepare one remote 1000G request per chromosome for all needed positions.
ref_base <- paste0(
  "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/",
  "1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
)
ref_requests <- lapply(sort(unique(pvar_target$CHROM)), function(chrom) {
  chrom_num <- sub("^chr", "", chrom)
  pv <- pvar_target[CHROM == chrom]
  bed <- file.path(ref_dir, sprintf("%s_target_positions.bed", chrom))
  fwrite(unique(pv[, .(CHROM, start0 = POS - 1L, end1 = POS)]), bed,
         sep = "\t", col.names = FALSE)
  data.table(
    chr = chrom, n_positions = uniqueN(pv$POS), bed_file = bed,
    source_url = sprintf("%s/1kGP_high_coverage_Illumina.chr%s.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
                         ref_base, chrom_num),
    bcf_file = file.path(ref_dir, sprintf("%s.1000G_highcov.EUR_unrelated.bcf", chrom)),
    gt_file = file.path(ref_dir, sprintf("%s.1000G_highcov.EUR_unrelated.GT.tsv.gz", chrom))
  )
})
ref_manifest <- rbindlist(ref_requests)

fwrite(targets, file.path(out_dir, "target_manifest.tsv"), sep = "\t")
fwrite(gene_context, file.path(out_dir, "gene_context_manifest.tsv"), sep = "\t")
fwrite(gene_intervals, file.path(out_dir, "gene_interval_manifest.tsv"), sep = "\t")
fwrite(regions, file.path(out_dir, "region_manifest.tsv"), sep = "\t")
fwrite(gwas_manifest, file.path(out_dir, "gwas_manifest.tsv"), sep = "\t")
fwrite(ref_manifest, file.path(out_dir, "reference_request_manifest.tsv"), sep = "\t")

summary <- data.table(
  target_rows = nrow(targets), mapk3_completed_rows = sum(targets$mapk3_pilot_complete),
  remaining_rows = sum(!targets$mapk3_pilot_complete), unique_genes = uniqueN(targets$gene_id),
  gene_context_fits = nrow(gene_context), collapsed_regions = nrow(regions),
  disorder_gene_pairs = uniqueN(targets[, paste(disorder, gene_id)])
)
fwrite(summary, file.path(out_dir, "manifest_summary.tsv"), sep = "\t")
print(summary)
