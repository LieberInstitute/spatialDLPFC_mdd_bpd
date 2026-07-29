#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(here)
})

here::i_am(".git/HEAD")

repo_root <- here()
code_dir <- here("code", "11_eQTL_coloc")
source(file.path(code_dir, "utils.R"), chdir = FALSE)

tables_dir <- here("processed-data", "11_eQTL_coloc", "seurat", "tables")
plink2_prefix <- here("processed-data", "00_genotypes", "plink2", "merged_maf05")

context_order <- names(SEURAT_CONTEXT_TO_DATASET_ID)
split_order <- "all"
gwasx_disorders <- DEFAULT_GWASX_DISORDERS

read_gzip_table <- function(path, ...) {
  if (!file.exists(path)) stop("Missing required file: ", path)
  fread(cmd = paste("gzip -dc", shQuote(path)), ...)
}

strict_unified <- read_gzip_table(file.path(tables_dir, "map_significant_unified.csv.gz"))
strict_pairs <- read_gzip_table(file.path(tables_dir, "map_significant_pairs.csv.gz"))

## load matched GWAS broadly enough for the p < 1e-5 exploratory cutoff.
gwas_matched <- load_matched_gwas_by_disorder(
  disorders = gwasx_disorders,
  plink2_prefix = plink2_prefix,
  repo_root = repo_root,
  si_min = GWAS_MATCH_SI_MIN
)
gwasx_stats <- gwasx_matched_from_gwas(gwas_matched)
degs <- load_DEGs(repo_root = repo_root, mode = "standard")

map_significant_unified_GWASx <- append_disorder_deg_flags(
  apply_gwasx_annotations(strict_unified, gwasx_stats, repo_root = repo_root),
  degs$disorder_related,
  disorders = gwasx_disorders
)
map_significant_pairs_GWASx <- append_disorder_deg_flags(
  apply_gwasx_annotations(strict_pairs, gwasx_stats, repo_root = repo_root),
  degs$disorder_related,
  disorders = gwasx_disorders
)
map_significant_summary_GWASx <- summarize_gwasx_significant_pairs(
  map_significant_pairs_GWASx,
  context_order = context_order,
  split_order = split_order,
  disorder_related = degs$disorder_related,
  gwasx_disorders = gwasx_disorders
)

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
out_unified <- file.path(tables_dir, "map_significant_unified_GWASx.csv.gz")
out_pairs <- file.path(tables_dir, "map_significant_pairs_GWASx.csv.gz")
out_summary <- file.path(tables_dir, "map_significant_summary_GWASx.csv")

fwrite(map_significant_unified_GWASx, out_unified)
fwrite(map_significant_pairs_GWASx, out_pairs)
fwrite(map_significant_summary_GWASx, out_summary)

## retain legacy aliases for existing downstream consumers.
legacy_unified <- file.path(tables_dir, "map_sig_unified_GWAS_relaxed.csv.gz")
legacy_pairs <- file.path(tables_dir, "map_sig_pairs_GWAS_relaxed.csv.gz")
legacy_summary <- file.path(tables_dir, "map_sig_summary_GWAS_relaxed.csv")
fwrite(map_significant_unified_GWASx, legacy_unified)
fwrite(map_significant_pairs_GWASx, legacy_pairs)
fwrite(map_significant_summary_GWASx, legacy_summary)

audit <- list(
  gwas_releases = GWAS_RELEASE_TAGS[gwasx_disorders],
  gwas_si_min = GWAS_MATCH_SI_MIN[gwasx_disorders],
  gwasx_thresholds = GWAS_EXPLORATORY_P_THRESHOLDS[gwasx_disorders],
  gwasx_variant_counts = vapply(gwasx_stats, nrow, integer(1)),
  outputs = c(
    map_significant_unified_GWASx = out_unified,
    map_significant_pairs_GWASx = out_pairs,
    map_significant_summary_GWASx = out_summary,
    map_sig_unified_GWAS_relaxed = legacy_unified,
    map_sig_pairs_GWAS_relaxed = legacy_pairs,
    map_sig_summary_GWAS_relaxed = legacy_summary
  ),
  broad_deg_gwasx_exact_genes = lapply(gwasx_disorders, function(dis) {
    flag_col <- paste0(dis, "_GWASx")
    sort(unique(map_significant_pairs_GWASx[DEG == 1L & get(flag_col) == 1L, gene_name]))
  })
)
names(audit$broad_deg_gwasx_exact_genes) <- gwasx_disorders

print(audit)
