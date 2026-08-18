#!/usr/bin/env Rscript

get_script_path <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", ca, value = TRUE)
  if (length(m) > 0) return(sub("^--file=", "", m[[1]]))

  for (i in rev(seq_len(sys.nframe()))) {
    of <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(of)) return(of)
  }

  stop("Cannot determine script path")
}

JHPCE_HOST <- "jh"
JHPCE_REPO_ROOT <- "/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd"

DEG_FILE_SPECS <- list(
  smoothed_layer_adjusted = file.path(
    "processed-data", "07_dx_DE",
    "layer-adjusted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv"
  ),
  smoothed_layer_restricted = file.path(
    "processed-data", "07_dx_DE",
    "layer-restricted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv"
  ),
  seurat_layer_adjusted = file.path(
    "processed-data", "07_dx_DE",
    "layer-adjusted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv"
  ),
  seurat_layer_restricted = file.path(
    "processed-data", "07_dx_DE",
    "layer-restricted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv"
  )
)

CUSTOM_CLUSTER_DEG_FILE_SPECS <- list(
  custom_layer_adjusted = file.path(
    "processed-data", "07_dx_DE",
    "layer-adjusted-pc3-age-nspots_custom-cluster_dx-sex_degs-F-test-t-test.csv"
  ),
  custom_layer_restricted = file.path(
    "processed-data", "07_dx_DE",
    "layer-restricted-pc3-age-nspots_custom-cluster_dx-sex_degs-F-test-t-test.csv"
  )
)

AUTHOR_UNION_REL_PATH <- file.path(
  "raw-data", "SCENIC_aux", "tf_lists", "MBv_PRECAST-Seurat_F-test-adjp-05.txt"
)

GENE_RANGES_REL_PATH <- file.path("processed-data", "ref", "granges.qs2")

GWAS_BCF_FILES <- c(
  BD = file.path("BD", "bip2024_eur_no23andMe.hg38.bcf"),
  MDD = file.path("MDD", "pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.bcf"),
  SCZD = file.path("SCZD", "PGC3_SCZ_wave3.european.autosome.public.v3.hg38.bcf")
)

DEFAULT_GWAS_OVERLAP_DISORDERS <- c("SCZD", "MDD", "BD")

GWAS_STRICT_P_THRESHOLD <- 5e-8

GWAS_MATCH_SI_MIN <- 0.8

## GWASx is the suggestive/exploratory mood-disorder overlap threshold.
## It is relaxed relative to strict genome-wide significance.
GWAS_EXPLORATORY_P_THRESHOLDS <- c(MDD = 1e-5, BD = 1e-5)

DEFAULT_GWASX_DISORDERS <- names(GWAS_EXPLORATORY_P_THRESHOLDS)

GWAS_GENE_LIST_FILES <- list(
  BD = list(
    broad = file.path("BD", "bd2024_gene_lists.tsv"),
    prio = file.path("BD", "bd2024_prioritized_credible_genes.tsv")
  ),
  MDD = list(
    broad = file.path("MDD", "mdd2025_high_confidence_genes.tsv"),
    prio = file.path("MDD", "mdd2025_high_confidence_genes.tsv")
  ),
  SCZD = list(
    broad = file.path("SCZD", "sczd2022_gene_lists.tsv"),
    prio = file.path("SCZD", "sczd2022_prioritized_genes.tsv")
  )
)

GWAS_STANDARD_GENE_LIST_FILES <- list(
  BD = list(
    broad = file.path("BD", "GWAS_BD_gene_list.tsv"),
    prio = file.path("BD", "GWAS_BD_prio_gene_list.tsv")
  ),
  MDD = list(
    broad = file.path("MDD", "GWAS_MDD_gene_list.tsv"),
    prio = file.path("MDD", "GWAS_MDD_prio_gene_list.tsv")
  ),
  SCZD = list(
    broad = file.path("SCZD", "GWAS_SCZD_gene_list.tsv"),
    prio = file.path("SCZD", "GWAS_SCZD_prio_gene_list.tsv")
  )
)

DEFAULT_GENE_RANGES_SOURCE_RELS <- c(
  file.path(
    "processed-data", "06_pseudobulk", "Seurat",
    "spe_n119_pseudo-no-lowUMI_sample-seurat-pc30_norm-filt.Rdata"
  )
)

SEURAT_CONTEXT_TO_DATASET_ID <- c(
  "Astro" = "astro",
  "Inhb" = "inhb",
  "L2.3" = "l2-3",
  "L4" = "l4",
  "L5" = "l5",
  "L6" = "l6",
  "Micro.Vasc" = "uvasc",
  "Oligo" = "oligo"
)

CUSTOM_CONTEXT_TO_DATASET_ID <- c(
  "Astro.L1" = "astro-l1",
  "Astro.Nrn" = "astro-nrn",
  "Inhb" = "inhb",
  "L2" = "l2",
  "L3" = "l3",
  "L4" = "l4",
  "L5" = "l5",
  "L6" = "l6",
  "Micro.Vasc" = "uvasc",
  "WM" = "wm"
)

DEG_SEX_PREFIX <- c(female = "F", male = "M")

## DEG t-test contrast tokens must match the 07_dx_DE summary column names, which
## use the legacy "BPD" label. Only the disorder key is the project-standard "BD".
DEG_TTEST_CONTRASTS <- c("NTC.MDD", "NTC.BPD", "MDD.BPD")

DEG_DISORDER_CONTRASTS <- c(MDD = "NTC.MDD", BD = "NTC.BPD")

DEG_TTEST_SIG_LABELS <- c("padj<.0001", "padj<.01", "padj<.05")

DEG_SPLITS <- c("all", "male", "female")

get_repo_root <- function(script_dir, args) {
  if (length(args) > 0) {
    normalizePath(args[[1]], mustWork = TRUE)
  } else {
    normalizePath(file.path(script_dir, "../.."), mustWork = TRUE)
  }
}

resolve_repo_root <- function(repo_root = NULL) {
  if (!is.null(repo_root)) return(normalizePath(repo_root, mustWork = TRUE))

  script_path <- tryCatch(get_script_path(), error = function(e) NULL)
  if (!is.null(script_path)) {
    script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
    return(get_repo_root(script_dir, commandArgs(trailingOnly = TRUE)))
  }

  if (file.exists(file.path(getwd(), "utils.R"))) {
    ## support source("utils.R") from this project directory.
    return(normalizePath(file.path(getwd(), "../.."), mustWork = TRUE))
  }

  stop("Cannot determine repo root; pass repo_root")
}

gwas_check_disorder <- function(dis) {
  dis <- toupper(dis)
  if (identical(dis, "BPD")) dis <- "BD"  # accept legacy "BPD" alias on input
  if (!dis %in% names(GWAS_BCF_FILES)) {
    stop("Unsupported disorder: ", dis, ". Expected one of: ", paste(names(GWAS_BCF_FILES), collapse = ", "))
  }
  dis
}

gwas_pval_tag <- function(pval) {
  if (length(pval) != 1 || is.na(pval) || !is.finite(pval) || pval <= 0 || pval >= 1) {
    stop("pval must be one finite value between 0 and 1")
  }
  tag <- formatC(pval, format = "e", digits = 0)
  tag <- sub("e-0+", "e-", tag)
  tag <- sub("e\\+0+", "e", tag)
  tag <- sub("e\\+", "e", tag)
  tag
}

gwas_si_tag <- function(si_min) {
  if (length(si_min) != 1 || is.na(si_min) || !is.finite(si_min)) {
    stop("si_min must be one finite value")
  }
  format(si_min, scientific = FALSE, trim = TRUE)
}

gwas_genotype_dir <- function(repo_root = NULL, genotype_dir = NULL) {
  if (!is.null(genotype_dir)) return(normalizePath(genotype_dir, mustWork = TRUE))
  file.path(resolve_repo_root(repo_root), "processed-data", "ref", "GWAS")
}

gwas_bcf_path <- function(dis, repo_root = NULL, genotype_dir = NULL) {
  dis <- gwas_check_disorder(dis)
  bcf_file <- file.path(gwas_genotype_dir(repo_root, genotype_dir), GWAS_BCF_FILES[[dis]])
  if (!file.exists(bcf_file)) stop("Missing GWAS BCF for ", dis, ": ", bcf_file)
  normalizePath(bcf_file, mustWork = TRUE)
}

gwas_cache_file <- function(dis, pval, repo_root = NULL, genotype_dir = NULL, si_min = 0.8) {
  dis <- gwas_check_disorder(dis)
  file.path(
    gwas_genotype_dir(repo_root, genotype_dir),
    sprintf("GWAS-%s_flt_p%s_SI%s.hg38.tab.gz", dis, gwas_pval_tag(pval), gwas_si_tag(si_min))
  )
}

gwas_gene_list_role <- function(use_prio = FALSE) {
  if (isTRUE(use_prio)) "prio" else "broad"
}

gwas_gene_list_rel <- function(dis, use_prio = FALSE, standardized = TRUE) {
  dis <- gwas_check_disorder(dis)
  role <- gwas_gene_list_role(use_prio)
  rels <- if (isTRUE(standardized)) GWAS_STANDARD_GENE_LIST_FILES else GWAS_GENE_LIST_FILES
  rel <- rels[[dis]][[role]]
  if (is.null(rel) || !nzchar(rel)) stop("Missing GWAS gene-list path for ", dis, " role ", role)
  rel
}

gwas_gene_list_path <- function(dis, use_prio = FALSE, repo_root = NULL,
                                genotype_dir = NULL, standardized = TRUE) {
  file.path(
    gwas_genotype_dir(repo_root = repo_root, genotype_dir = genotype_dir),
    gwas_gene_list_rel(dis, use_prio = use_prio, standardized = standardized)
  )
}

gwas_collapse_values <- function(x) {
  x <- sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))])))
  if (length(x) == 0) "" else paste(x, collapse = ";")
}

standardize_gwas_gene_list <- function(dt, dis, use_prio = FALSE) {
  dis <- gwas_check_disorder(dis)
  if (!"gene_symbol" %in% names(dt)) stop("Missing gene_symbol column for ", dis, " GWAS gene list")
  dt <- data.table::as.data.table(dt)
  dt <- dt[!is.na(gene_symbol) & nzchar(gene_symbol)]
  if (nrow(dt) == 0) {
    return(data.table::data.table(
      disorder = character(),
      gene_symbol = character(),
      list_role = character(),
      source_list_names = character(),
      primary_list = character(),
      evidence_types = character(),
      source_files = character(),
      n_source_rows = integer()
    ))
  }

  role <- if (isTRUE(use_prio)) "prio" else "broad"
  out <- dt[, .(
    disorder = dis,
    list_role = role,
    source_list_names = if ("list_name" %in% names(.SD)) gwas_collapse_values(list_name) else "",
    primary_list = if ("primary_list" %in% names(.SD) && any(primary_list == "yes", na.rm = TRUE)) "yes" else "no",
    evidence_types = if ("evidence_type" %in% names(.SD)) gwas_collapse_values(evidence_type) else "",
    source_files = if ("source_file" %in% names(.SD)) gwas_collapse_values(source_file) else "",
    n_source_rows = .N
  ), by = .(gene_symbol)]
  data.table::setcolorder(
    out,
    c("disorder", "gene_symbol", "list_role", "source_list_names", "primary_list",
      "evidence_types", "source_files", "n_source_rows")
  )
  out[order(gene_symbol)]
}

write_standard_gwas_gene_lists <- function(repo_root = NULL, genotype_dir = NULL,
                                           disorders = DEFAULT_GWAS_OVERLAP_DISORDERS) {
  out <- list()
  for (dis in disorders) {
    dis <- gwas_check_disorder(dis)
    for (use_prio in c(FALSE, TRUE)) {
      src <- gwas_gene_list_path(dis, use_prio = use_prio, repo_root = repo_root,
                                 genotype_dir = genotype_dir, standardized = FALSE)
      dest <- gwas_gene_list_path(dis, use_prio = use_prio, repo_root = repo_root,
                                  genotype_dir = genotype_dir, standardized = TRUE)
      if (!file.exists(src)) stop("Missing source GWAS gene list: ", src)
      dt <- data.table::fread(src)
      std <- standardize_gwas_gene_list(dt, dis = dis, use_prio = use_prio)
      dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
      data.table::fwrite(std, dest, sep = "\t", quote = FALSE)
      out[[paste(dis, gwas_gene_list_role(use_prio), sep = "_")]] <- dest
    }
  }
  unlist(out, use.names = TRUE)
}

loadGWASGeneList <- function(dis, use_prio = FALSE, repo_root = NULL,
                             genotype_dir = NULL, prefer_standard = TRUE) {
  dis <- gwas_check_disorder(dis)
  std_path <- gwas_gene_list_path(dis, use_prio = use_prio, repo_root = repo_root,
                                  genotype_dir = genotype_dir, standardized = TRUE)
  if (isTRUE(prefer_standard) && file.exists(std_path)) {
    dt <- data.table::fread(std_path)
  } else {
    src_path <- gwas_gene_list_path(dis, use_prio = use_prio, repo_root = repo_root,
                                    genotype_dir = genotype_dir, standardized = FALSE)
    if (!file.exists(src_path)) stop("Missing source GWAS gene list: ", src_path)
    dt <- standardize_gwas_gene_list(data.table::fread(src_path), dis = dis, use_prio = use_prio)
  }

  req <- c("disorder", "gene_symbol", "list_role", "source_list_names", "primary_list",
           "evidence_types", "source_files", "n_source_rows")
  missing_cols <- setdiff(req, names(dt))
  if (length(missing_cols) > 0) {
    stop("Missing standardized GWAS gene-list columns: ", paste(missing_cols, collapse = ", "))
  }
  dt[!is.na(gene_symbol) & nzchar(gene_symbol)]
}

normalize_gwas_query_table <- function(dt, dis, pval, cache_file = NULL, si_min = 0.8) {
  dis <- gwas_check_disorder(dis)
  query_cols <- c("chr", "pos", "rsid", "a0", "a1", "beta", "beta_se", "lp", "N", "ns", "ncas", "impinfo")

  if (nrow(dt) == 0) {
    dt <- data.table::as.data.table(stats::setNames(replicate(length(query_cols), logical(), simplify = FALSE), query_cols))
  } else {
    if (ncol(dt) != length(query_cols)) {
      stop("Unexpected GWAS query column count: ", ncol(dt), "; expected ", length(query_cols))
    }
    data.table::setnames(dt, query_cols)
  }

  dt[, `:=`(
    chr = as.character(chr),
    pos = as.integer(pos),
    rsid = as.character(rsid),
    a0 = as.character(a0),
    a1 = as.character(a1),
    beta = as.numeric(beta),
    beta_se = as.numeric(beta_se),
    lp = as.numeric(lp),
    N = as.numeric(N),
    ns = as.numeric(ns),
    ncas = as.numeric(ncas),
    impinfo = as.numeric(impinfo)
  )]
  dt[, p := 10^(-lp)]
  dt[, variant_id := sprintf("%s:%s:%s:%s", chr, pos, a0, a1)]
  data.table::setcolorder(dt, c("rsid", "chr", "pos", "a0", "a1", "beta", "beta_se", "N", "p", "impinfo", "ncas", "ns", "lp", "variant_id"))

  attr(dt, "gwas_dis") <- dis
  attr(dt, "gwas_pval") <- pval
  attr(dt, "gwas_si_min") <- si_min
  if (!is.null(cache_file)) {
    attr(dt, "gwas_cache_file") <- cache_file
    attr(dt, "gwas_cache_dir") <- dirname(cache_file)
  }
  dt
}

annotate_gwas_table <- function(dt, dis, pval, cache_file, si_min = 0.8) {
  attr(dt, "gwas_dis") <- gwas_check_disorder(dis)
  attr(dt, "gwas_pval") <- pval
  attr(dt, "gwas_si_min") <- si_min
  attr(dt, "gwas_cache_file") <- cache_file
  attr(dt, "gwas_cache_dir") <- dirname(cache_file)
  dt
}

loadGWAS <- function(dis, pval, repo_root = NULL, genotype_dir = NULL,
                     si_min = 0.8, bcftools = "bcftools", use_cache = TRUE) {
  dis <- gwas_check_disorder(dis)
  cache_file <- gwas_cache_file(dis, pval, repo_root = repo_root, genotype_dir = genotype_dir, si_min = si_min)
  if (use_cache && file.exists(cache_file)) {
    return(annotate_gwas_table(data.table::fread(cache_file), dis, pval, cache_file, si_min))
  }

  bcf_file <- gwas_bcf_path(dis, repo_root = repo_root, genotype_dir = genotype_dir)
  lp_min <- -log10(pval)
  query_file <- tempfile(pattern = "gwas-query-", fileext = ".tab")
  err_file <- tempfile(pattern = "gwas-query-err-", fileext = ".log")
  on.exit(unlink(c(query_file, err_file)), add = TRUE)

  ## query only the fields required by downstream eQTL and coloc steps.
  status <- system2(
    bcftools,
    args = c(
      "query",
      "-i", shQuote(sprintf("FORMAT/LP>=%s && FORMAT/SI>=%s", format(lp_min, scientific = FALSE), gwas_si_tag(si_min))),
      "-f", shQuote("%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%ES\t%SE\t%LP\t%NE\t%NS\t%NC\t%SI]\n"),
      shQuote(bcf_file)
    ),
    stdout = query_file,
    stderr = err_file
  )
  if (!identical(status, 0L)) {
    err <- if (file.exists(err_file)) readLines(err_file, warn = FALSE) else character()
    stop(
      "bcftools query failed for ", dis, ": ", bcf_file,
      if (length(err) > 0) paste0("\n", paste(err, collapse = "\n")) else ""
    )
  }

  dt <- if (file.exists(query_file) && file.info(query_file)$size > 0) {
    data.table::fread(query_file, header = FALSE)
  } else {
    data.table::data.table()
  }
  dt <- normalize_gwas_query_table(dt, dis = dis, pval = pval, cache_file = cache_file, si_min = si_min)

  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(dt, cache_file, sep = "\t", quote = FALSE)
  annotate_gwas_table(dt, dis, pval, cache_file, si_min)
}

gwas_tqtl_cache_file <- function(gwas, plink2_prefix) {
  dis <- attr(gwas, "gwas_dis")
  pval <- attr(gwas, "gwas_pval")
  si_min <- attr(gwas, "gwas_si_min")
  cache_dir <- attr(gwas, "gwas_cache_dir")
  if (is.null(dis) || is.null(pval) || is.null(si_min) || is.null(cache_dir)) {
    stop("gwas must be returned by loadGWAS() so cache metadata is available")
  }
  file.path(
    cache_dir,
    sprintf(
      "GWAS-%s_flt_p%s_SI%s_%s_tqtl-matched.tab.gz",
      dis, gwas_pval_tag(pval), gwas_si_tag(si_min), basename(plink2_prefix)
    )
  )
}

read_plink2_pvar <- function(plink2_prefix) {
  pvar_file <- paste0(plink2_prefix, ".pvar")
  if (!file.exists(pvar_file)) stop("Missing PLINK2 pvar file: ", pvar_file)
  pvar <- data.table::fread(pvar_file, skip = "#CHROM")
  if ("#CHROM" %in% names(pvar)) data.table::setnames(pvar, "#CHROM", "CHROM")
  req <- c("CHROM", "POS", "ID", "REF", "ALT")
  missing_cols <- setdiff(req, names(pvar))
  if (length(missing_cols) > 0) {
    stop("Missing required pvar columns: ", paste(missing_cols, collapse = ", "))
  }
  pvar[, .(
    chr = as.character(CHROM),
    pos = as.integer(POS),
    pvar_ref = as.character(REF),
    pvar_alt = as.character(ALT),
    pvar_variant_id = as.character(ID)
  )]
}

genotype_rsid_vcf_path <- function(repo_root = NULL) {
  file.path(
    resolve_repo_root(repo_root),
    "processed-data", "00_genotypes", "merged_R.8_MAF.01.RSann.vcf.gz"
  )
}

genotype_variant_info_cache_file <- function(plink2_prefix) {
  paste0(plink2_prefix, "_variant_info.csv.gz")
}

read_plink2_variant_table <- function(plink2_prefix) {
  pvar_file <- paste0(plink2_prefix, ".pvar")
  if (!file.exists(pvar_file)) stop("Missing PLINK2 pvar file: ", pvar_file)
  pvar <- data.table::fread(pvar_file, skip = "#CHROM")
  if ("#CHROM" %in% names(pvar)) data.table::setnames(pvar, "#CHROM", "CHROM")
  req <- c("CHROM", "POS", "ID", "REF", "ALT")
  missing_cols <- setdiff(req, names(pvar))
  if (length(missing_cols) > 0) {
    stop("Missing required pvar columns: ", paste(missing_cols, collapse = ", "))
  }

  out <- pvar[, .(
    variant_id = as.character(ID),
    CHROM = as.character(CHROM),
    POS = as.integer(POS),
    REF = as.character(REF),
    ALT = as.character(ALT),
    var_idx = .I
  )]
  if (anyDuplicated(out$variant_id)) stop("Duplicate pvar variant IDs in ", pvar_file)
  out
}

query_genotype_vcf_rsids <- function(pvar_info, rsid_vcf, bcftools = "bcftools") {
  if (!nzchar(Sys.which(bcftools))) stop("Missing bcftools executable: ", bcftools)
  if (!file.exists(rsid_vcf)) stop("Missing genotype rsID VCF: ", rsid_vcf)

  region_file <- tempfile("genotype-rsid-regions-", fileext = ".bed")
  query_file <- tempfile("genotype-rsid-query-", fileext = ".tsv")
  err_file <- tempfile("genotype-rsid-query-err-", fileext = ".log")
  on.exit(unlink(c(region_file, query_file, err_file)), add = TRUE)

  ## query only pvar positions; exact REF/ALT matching happens after query.
  regions <- unique(pvar_info[, .(CHROM, start0 = POS - 1L, end1 = POS)])
  data.table::fwrite(regions, region_file, sep = "\t", col.names = FALSE)

  status <- system2(
    bcftools,
    args = c(
      "query",
      "-R", shQuote(region_file),
      "-f", shQuote("%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/RS\\n"),
      shQuote(rsid_vcf)
    ),
    stdout = query_file,
    stderr = err_file
  )
  if (!identical(status, 0L)) {
    err <- if (file.exists(err_file)) readLines(err_file, warn = FALSE) else character()
    stop(
      "bcftools rsID query failed for genotype VCF: ", rsid_vcf,
      if (length(err) > 0) paste0("\n", paste(err, collapse = "\n")) else ""
    )
  }

  if (!file.exists(query_file) || file.info(query_file)$size == 0) {
    return(data.table::data.table(
      variant_id = character(),
      rs_numeric = character(),
      rsid = character(),
      rsid_source = character()
    ))
  }

  q <- data.table::fread(
    query_file,
    header = FALSE,
    col.names = c("CHROM", "POS", "vcf_id", "REF", "ALT", "RS")
  )
  q[, `:=`(
    variant_id = paste(CHROM, POS, REF, ALT, sep = ":"),
    rs_numeric = data.table::fifelse(is.na(RS) | RS == "." | RS == "", "", as.character(RS))
  )]
  q[, rsid := data.table::fifelse(nzchar(rs_numeric), paste0("rs", rs_numeric), "")]
  q[, rsid_source := data.table::fifelse(nzchar(rsid), "genotype_vcf_INFO_RS", "")]
  q[!is.na(variant_id) & nzchar(variant_id), .(
    rs_numeric = rs_numeric[match(TRUE, nzchar(rs_numeric), nomatch = 1L)],
    rsid = rsid[match(TRUE, nzchar(rsid), nomatch = 1L)],
    rsid_source = rsid_source[match(TRUE, nzchar(rsid), nomatch = 1L)]
  ), by = variant_id]
}

build_genotype_variant_info_cache <- function(plink2_prefix,
                                              rsid_vcf = genotype_rsid_vcf_path(),
                                              out_file = genotype_variant_info_cache_file(plink2_prefix),
                                              bcftools = "bcftools") {
  pvar_info <- read_plink2_variant_table(plink2_prefix)
  rsids <- query_genotype_vcf_rsids(pvar_info, rsid_vcf = rsid_vcf, bcftools = bcftools)
  data.table::setkey(rsids, variant_id)
  out <- rsids[pvar_info, on = "variant_id"]
  out[is.na(rsid), `:=`(rsid = "", rs_numeric = "", rsid_source = "")]
  out <- out[, .(variant_id, CHROM, POS, REF, ALT, var_idx, rsid, rs_numeric, rsid_source)]
  data.table::setorder(out, var_idx)

  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(out, out_file)
  out
}

load_genotype_variant_info <- function(plink2_prefix,
                                       rsid_vcf = genotype_rsid_vcf_path(),
                                       cache_file = genotype_variant_info_cache_file(plink2_prefix),
                                       use_cache = TRUE,
                                       build_if_missing = TRUE,
                                       bcftools = "bcftools") {
  if (isTRUE(use_cache) && file.exists(cache_file)) {
    out <- data.table::fread(cache_file)
  } else {
    if (!isTRUE(build_if_missing)) {
      stop(
        "Missing genotype variant cache: ", cache_file,
        "\nBuild it in 03b_eQTL_boxplots.Rmd or copy it from the machine where it was built."
      )
    }
    out <- build_genotype_variant_info_cache(
      plink2_prefix = plink2_prefix,
      rsid_vcf = rsid_vcf,
      out_file = cache_file,
      bcftools = bcftools
    )
  }

  req <- c("variant_id", "CHROM", "POS", "REF", "ALT", "var_idx", "rsid", "rs_numeric", "rsid_source")
  missing_cols <- setdiff(req, names(out))
  if (length(missing_cols) > 0) {
    stop("Missing genotype variant cache columns: ", paste(missing_cols, collapse = ", "))
  }
  if (anyDuplicated(out$variant_id)) stop("Duplicate variant IDs in genotype variant cache: ", cache_file)
  out[, `:=`(
    rsid = as.character(rsid),
    rs_numeric = as.character(rs_numeric),
    rsid_source = as.character(rsid_source)
  )]
  out[is.na(rsid), rsid := ""]
  out[is.na(rs_numeric), rs_numeric := ""]
  out[is.na(rsid_source), rsid_source := ""]
  out
}

matchGwasGeno <- function(gwas, plink2_prefix, use_cache = TRUE) {
  cache_file <- gwas_tqtl_cache_file(gwas, plink2_prefix)
  if (use_cache && file.exists(cache_file)) {
    return(data.table::fread(cache_file))
  }

  req <- c("rsid", "chr", "pos", "a0", "a1", "beta", "beta_se", "N", "p", "ncas", "impinfo")
  missing_cols <- setdiff(req, names(gwas))
  if (length(missing_cols) > 0) {
    stop("Missing required GWAS columns: ", paste(missing_cols, collapse = ", "))
  }

  pvar_info <- read_plink2_pvar(plink2_prefix)
  gwas_dt <- data.table::copy(gwas)
  gwas_dt[, gwas_row_id := .I]

  exact <- merge(
    gwas_dt,
    pvar_info,
    by.x = c("chr", "pos", "a0", "a1"),
    by.y = c("chr", "pos", "pvar_ref", "pvar_alt"),
    all = FALSE,
    allow.cartesian = TRUE
  )
  if (nrow(exact) > 0) {
    exact[, `:=`(
      variant_id = pvar_variant_id,
      A1 = a1,
      A2 = a0,
      match_mode = "exact",
      match_rank = 1L
    )]
  }

  swapped <- merge(
    gwas_dt,
    pvar_info,
    by.x = c("chr", "pos", "a0", "a1"),
    by.y = c("chr", "pos", "pvar_alt", "pvar_ref"),
    all = FALSE,
    allow.cartesian = TRUE
  )
  if (nrow(swapped) > 0) {
    swapped[, `:=`(
      beta = -beta,
      variant_id = pvar_variant_id,
      A1 = a0,
      A2 = a1,
      match_mode = "swapped",
      match_rank = 2L
    )]
  }

  gwas_tqtl <- data.table::rbindlist(list(exact, swapped), use.names = TRUE, fill = TRUE)
  out_cols <- c("variant_id", "rsid", "A1", "A2", "beta", "beta_se", "N", "p", "ncas", "impinfo", "match_mode")
  if (nrow(gwas_tqtl) == 0) {
    gwas_tqtl <- data.table::as.data.table(stats::setNames(replicate(length(out_cols), logical(), simplify = FALSE), out_cols))
  } else {
    ## prefer exact matches if duplicate variant IDs appear.
    data.table::setorder(gwas_tqtl, variant_id, match_rank, p)
    gwas_tqtl <- gwas_tqtl[!duplicated(variant_id), ..out_cols]
  }

  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(gwas_tqtl, cache_file, sep = "\t", quote = FALSE)
  gwas_tqtl
}

fetch_remote_file <- function(host, remote_path, local_path) {
  dir.create(dirname(local_path), recursive = TRUE, showWarnings = FALSE)

  test_status <- system2(
    "ssh",
    args = c("-o", "BatchMode=yes", host, "test", "-r", remote_path),
    stdout = FALSE,
    stderr = FALSE
  )
  if (!identical(test_status, 0L)) {
    stop("Remote DEG file is not readable via ", host, ": ", remote_path)
  }

  tmp_file <- tempfile(pattern = "deg-fetch-", tmpdir = dirname(local_path))
  err_file <- tempfile(pattern = "deg-fetch-err-", tmpdir = dirname(local_path))
  on.exit(unlink(c(tmp_file, err_file)), add = TRUE)

  status <- system2(
    "ssh",
    args = c("-o", "BatchMode=yes", host, "cat", remote_path),
    stdout = tmp_file,
    stderr = err_file
  )
  if (!identical(status, 0L)) {
    err <- character()
    if (file.exists(err_file)) {
      err <- readLines(err_file, warn = FALSE)
    }
    stop(
      "Failed to fetch DEG file from ", host, ": ", remote_path,
      if (length(err) > 0) paste0("\n", paste(err, collapse = "\n")) else ""
    )
  }
  if (!file.exists(tmp_file) || file.info(tmp_file)$size == 0) {
    stop("Fetched DEG file is empty: ", remote_path)
  }

  if (!file.rename(tmp_file, local_path)) {
    ok <- file.copy(tmp_file, local_path, overwrite = TRUE)
    if (!ok) stop("Failed to stage DEG file locally: ", local_path)
    unlink(tmp_file)
  }
}

ensure_local_file <- function(repo_root, rel_path, host = JHPCE_HOST, remote_root = JHPCE_REPO_ROOT) {
  local_path <- file.path(repo_root, rel_path)
  if (file.exists(local_path)) return(local_path)

  remote_path <- file.path(remote_root, rel_path)
  message("Staging missing file from ", host, ": ", rel_path)
  fetch_remote_file(host = host, remote_path = remote_path, local_path = local_path)
  local_path
}

resolve_gene_range_source_file <- function(path, repo_root, host, remote_root) {
  if (grepl("^/", path)) {
    if (!file.exists(path)) stop("Missing gene range source file: ", path)
    return(normalizePath(path, mustWork = TRUE))
  }

  ensure_local_file(
    repo_root = repo_root,
    rel_path = path,
    host = host,
    remote_root = remote_root
  )
}

load_summarized_experiment_objects <- function(path) {
  obj_env <- new.env(parent = emptyenv())
  obj_names <- load(path, envir = obj_env)
  if (length(obj_names) == 0) stop("No objects found in ", path)

  is_se <- vapply(
    obj_names,
    function(x) inherits(get(x, envir = obj_env), "SummarizedExperiment"),
    logical(1)
  )
  if (!any(is_se)) {
    stop("No SummarizedExperiment object found in ", path)
  }

  lapply(obj_names[is_se], function(x) get(x, envir = obj_env))
}

validate_gene_ranges <- function(granges, label) {
  req <- c("gene_id", "gene_name")
  missing_cols <- setdiff(req, names(S4Vectors::mcols(granges)))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required gene range metadata in ", label, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  if (is.null(names(granges)) || anyNA(names(granges)) || any(!nzchar(names(granges)))) {
    stop("Missing GRanges names in ", label)
  }

  invisible(TRUE)
}

load_gene_ranges <- function(repo_root = NULL, source_files = NULL, granges_file = NULL,
                             verbose = TRUE, host = JHPCE_HOST,
                             remote_root = JHPCE_REPO_ROOT) {
  repo_root <- resolve_repo_root(repo_root)

  if (is.null(granges_file)) {
    granges_file <- file.path(repo_root, GENE_RANGES_REL_PATH)
  }
  granges_file <- normalizePath(granges_file, mustWork = FALSE)

  if (file.exists(granges_file)) {
    granges <- qs2::qs_read(granges_file)
    validate_gene_ranges(granges, granges_file)
    if (isTRUE(verbose)) cat("Loaded gene ranges:", granges_file, "\n")
    return(granges)
  }

  if (is.null(source_files)) source_files <- DEFAULT_GENE_RANGES_SOURCE_RELS
  source_paths <- vapply(
    source_files,
    resolve_gene_range_source_file,
    character(1),
    repo_root = repo_root,
    host = host,
    remote_root = remote_root
  )

  se_objects <- unlist(lapply(source_paths, load_summarized_experiment_objects), recursive = FALSE)
  all_granges <- lapply(se_objects, SummarizedExperiment::rowRanges)
  merged_granges <- do.call(c, unname(all_granges))

  ## preserve first occurrence when multiple source objects share gene IDs.
  merged_granges <- merged_granges[!duplicated(names(merged_granges))]
  validate_gene_ranges(merged_granges, paste(source_paths, collapse = ", "))

  dir.create(dirname(granges_file), recursive = TRUE, showWarnings = FALSE)
  qs2::qs_save(merged_granges, file = granges_file)
  if (isTRUE(verbose)) cat("Built gene ranges:", granges_file, "\n")

  merged_granges
}

load_deg_summary <- function(path, label) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  req <- c("gene_id", "gene_name", "adj.P.Val")
  missing_cols <- setdiff(req, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required DEG columns in ", label, " (", path, "): ",
      paste(missing_cols, collapse = ", ")
    )
  }

  df[df$adj.P.Val < 0.05, , drop = FALSE]
}

read_author_union_list <- function(path) {
  genes <- scan(path, what = "character", quiet = TRUE)
  sort(unique(genes[nzchar(genes)]))
}

collapse_gene_table <- function(df, label) {
  req <- c("gene_id", "gene_name")
  missing_cols <- setdiff(req, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing gene key columns in ", label, ": ", paste(missing_cols, collapse = ", "))
  }

  df2 <- df[, req, drop = FALSE]
  split_names <- split(df2$gene_name, df2$gene_id)
  bad_ids <- names(split_names)[vapply(split_names, function(x) length(unique(x)) != 1L, logical(1))]
  if (length(bad_ids) > 0) {
    stop(
      "Inconsistent gene_name values for gene_id(s) in ", label, ": ",
      paste(head(bad_ids, 10), collapse = ", ")
    )
  }

  df2[!duplicated(df2$gene_id), , drop = FALSE]
}

dataset_id_map_from_context_sets <- function(context_sets, context_to_dataset_id = SEURAT_CONTEXT_TO_DATASET_ID) {
  out <- list()
  for (ctx in names(context_to_dataset_id)) {
    ds_id <- unname(context_to_dataset_id[[ctx]])
    set_df <- context_sets[[ctx]]
    out[[ds_id]] <- set_df
    out[[paste0(ds_id, "_m")]] <- set_df
    out[[paste0(ds_id, "_f")]] <- set_df
  }
  out
}

validate_author_union <- function(deg_global, author_union_genes) {
  global_gene_names <- sort(unique(deg_global$gene_name))
  missing_from_author_union <- setdiff(global_gene_names, author_union_genes)
  extra_in_author_union <- setdiff(author_union_genes, global_gene_names)

  list(
    matches = length(missing_from_author_union) == 0L && length(extra_in_author_union) == 0L,
    gene_name_count = length(global_gene_names),
    author_gene_name_count = length(author_union_genes),
    missing_from_author_union = missing_from_author_union,
    extra_in_author_union = extra_in_author_union
  )
}

ttest_sig_rows <- function(df, cols) {
  if (length(cols) == 0) return(rep(FALSE, nrow(df)))
  ## t-test summary files store adjusted p-value bins as strings.
  sig_df <- as.data.frame(
    lapply(df[, cols, drop = FALSE], function(x) x %in% DEG_TTEST_SIG_LABELS),
    check.names = FALSE
  )
  rowSums(sig_df, na.rm = TRUE) > 0
}

disorder_ttest_cols <- function(df, contrast) {
  contrast_regex <- gsub(".", "[.]", contrast, fixed = TRUE)
  grep(paste0("(^|_)[FM]_", contrast_regex, "_ttest$"), names(df), value = TRUE)
}

build_disorder_related_degs <- function(deg_global, deg_tables,
                                        disorder_contrasts = DEG_DISORDER_CONTRASTS) {
  broad_gene_ids <- unique(deg_global$gene_id)
  disorder_sets <- lapply(names(disorder_contrasts), function(disorder) {
    contrast <- unname(disorder_contrasts[[disorder]])
    out <- lapply(deg_tables, function(df) {
      cols <- disorder_ttest_cols(df, contrast)
      df[ttest_sig_rows(df, cols), c("gene_id", "gene_name"), drop = FALSE]
    })
    collapse_gene_table(do.call(rbind, out), label = paste0(disorder, "_related_DEGs"))
  })
  names(disorder_sets) <- names(disorder_contrasts)

  missing_from_broad <- lapply(disorder_sets, function(df) setdiff(df$gene_id, broad_gene_ids))
  list(
    global = disorder_sets,
    validation = list(
      global_gene_id_counts = vapply(disorder_sets, nrow, integer(1)),
      subset_of_broad = vapply(missing_from_broad, function(x) length(x) == 0L, logical(1)),
      missing_from_broad = missing_from_broad
    )
  )
}

validate_deg_view_columns <- function(seurat_lr_df,
                                      contexts = names(SEURAT_CONTEXT_TO_DATASET_ID),
                                      sex_prefix = DEG_SEX_PREFIX,
                                      contrasts = DEG_TTEST_CONTRASTS) {
  context_cols <- paste0("n_ttest_sig_", contexts)
  missing_context_cols <- setdiff(context_cols, names(seurat_lr_df))
  if (length(missing_context_cols) > 0) {
    stop("Missing Seurat context DEG columns: ", paste(missing_context_cols, collapse = ", "))
  }

  context_sex_cols <- unlist(lapply(contexts, function(context_name) {
    unlist(lapply(unname(sex_prefix), function(prefix) {
      paste0(context_name, "_", prefix, "_", contrasts, "_ttest")
    }), use.names = FALSE)
  }), use.names = FALSE)
  missing_context_sex_cols <- setdiff(context_sex_cols, names(seurat_lr_df))
  if (length(missing_context_sex_cols) > 0) {
    stop(
      "Missing Seurat context-by-sex DEG columns: ",
      paste(missing_context_sex_cols, collapse = ", ")
    )
  }

  invisible(TRUE)
}

make_deg_view_rows <- function(df, deg_view, context, split, deg_sex = "all") {
  df2 <- collapse_gene_table(df, label = paste(deg_view, context, split, deg_sex, sep = ":"))
  if (nrow(df2) == 0) {
    return(data.frame(
      deg_view = character(),
      context = character(),
      split = character(),
      deg_sex = character(),
      gene_id = character(),
      gene_name = character(),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    deg_view = deg_view,
    context = context,
    split = split,
    deg_sex = deg_sex,
    gene_id = df2$gene_id,
    gene_name = df2$gene_name,
    stringsAsFactors = FALSE
  )
}

sex_specific_deg_set <- function(deg_tables, sex, sex_prefix = DEG_SEX_PREFIX) {
  prefix <- unname(sex_prefix[[sex]])
  if (is.na(prefix)) stop("Unknown DEG sex: ", sex)

  out <- lapply(names(deg_tables), function(label) {
    df <- deg_tables[[label]]
    sex_cols <- grep(paste0("(^|_)", prefix, "_.*_ttest$"), names(df), value = TRUE)
    df[ttest_sig_rows(df, sex_cols), c("gene_id", "gene_name"), drop = FALSE]
  })

  collapse_gene_table(do.call(rbind, out), label = paste0("sex_specific:", sex))
}

build_deg_views <- function(deg_global, deg_tables,
                            contexts = names(SEURAT_CONTEXT_TO_DATASET_ID),
                            sex_prefix = DEG_SEX_PREFIX,
                            contrasts = DEG_TTEST_CONTRASTS,
                            splits = DEG_SPLITS) {
  seurat_lr_df <- deg_tables$seurat_layer_restricted
  validate_deg_view_columns(
    seurat_lr_df = seurat_lr_df,
    contexts = contexts,
    sex_prefix = sex_prefix,
    contrasts = contrasts
  )

  rows <- list()

  for (split_name in splits) {
    for (context_name in contexts) {
      rows[[length(rows) + 1L]] <- make_deg_view_rows(
        deg_global,
        deg_view = "broad_interaction",
        context = context_name,
        split = split_name
      )
    }
  }

  for (context_name in contexts) {
    col_name <- paste0("n_ttest_sig_", context_name)
    localized <- seurat_lr_df[
      !is.na(seurat_lr_df[[col_name]]) & seurat_lr_df[[col_name]] > 0,
      ,
      drop = FALSE
    ]
    for (split_name in splits) {
      rows[[length(rows) + 1L]] <- make_deg_view_rows(
        localized,
        deg_view = "context_localized",
        context = context_name,
        split = split_name
      )
    }
  }

  sex_sets <- lapply(names(sex_prefix), function(sex) {
    sex_specific_deg_set(deg_tables = deg_tables, sex = sex, sex_prefix = sex_prefix)
  })
  names(sex_sets) <- names(sex_prefix)

  for (context_name in contexts) {
    for (sex in names(sex_sets)) {
      for (split_name in c("all", sex)) {
        rows[[length(rows) + 1L]] <- make_deg_view_rows(
          sex_sets[[sex]],
          deg_view = "sex_specific",
          context = context_name,
          split = split_name,
          deg_sex = sex
        )
      }
    }
  }

  for (context_name in contexts) {
    for (sex in names(sex_prefix)) {
      prefix <- unname(sex_prefix[[sex]])
      sex_context_cols <- paste0(context_name, "_", prefix, "_", contrasts, "_ttest")
      sex_context <- seurat_lr_df[ttest_sig_rows(seurat_lr_df, sex_context_cols), , drop = FALSE]
      for (split_name in c("all", sex)) {
        rows[[length(rows) + 1L]] <- make_deg_view_rows(
          sex_context,
          deg_view = "context_and_sex_specific",
          context = context_name,
          split = split_name,
          deg_sex = sex
        )
      }
    }
  }

  long <- do.call(rbind, rows)
  rownames(long) <- NULL

  counts <- aggregate(
    gene_id ~ deg_view + context + split + deg_sex,
    data = long,
    FUN = function(x) length(unique(x)),
    na.action = NULL
  )
  names(counts)[names(counts) == "gene_id"] <- "n_genes"
  rownames(counts) <- NULL

  list(
    long = long,
    counts = counts,
    sex_sets = sex_sets
  )
}

build_custom_cluster_deg_views <- function(custom_by_context, splits = c("all")) {
  rows <- list()
  for (split_name in splits) {
    for (context_name in names(custom_by_context)) {
      rows[[length(rows) + 1L]] <- make_deg_view_rows(
        custom_by_context[[context_name]],
        deg_view = "custom_primary",
        context = context_name,
        split = split_name
      )
    }
  }

  long <- do.call(rbind, rows)
  rownames(long) <- NULL

  counts <- aggregate(
    gene_id ~ deg_view + context + split + deg_sex,
    data = long,
    FUN = function(x) length(unique(x)),
    na.action = NULL
  )
  names(counts)[names(counts) == "gene_id"] <- "n_genes"
  rownames(counts) <- NULL

  list(
    long = long,
    counts = counts,
    sex_sets = list()
  )
}

load_standard_DEGs <- function(repo_root = NULL, verbose = TRUE,
                      include_disorder_degs = TRUE,
                      host = JHPCE_HOST, remote_root = JHPCE_REPO_ROOT) {
  repo_root <- resolve_repo_root(repo_root)

  deg_file_paths <- lapply(DEG_FILE_SPECS, function(rel_path) {
    ensure_local_file(
      repo_root = repo_root,
      rel_path = rel_path,
      host = host,
      remote_root = remote_root
    )
  })
  author_union_file <- ensure_local_file(
    repo_root = repo_root,
    rel_path = AUTHOR_UNION_REL_PATH,
    host = host,
    remote_root = remote_root
  )

  deg_tables <- Map(
    f = function(path, label) load_deg_summary(path = path, label = label),
    path = deg_file_paths,
    label = names(deg_file_paths)
  )
  deg_tables <- deg_tables[names(DEG_FILE_SPECS)]

  ## global DEG support uses all four author-recommended F-test summaries.
  deg_global <- collapse_gene_table(
    do.call(
      rbind,
      lapply(unname(deg_tables), function(df) df[, c("gene_id", "gene_name"), drop = FALSE])
    ),
    label = "combined DEG summaries"
  )

  seurat_layer_adjusted_genes <- collapse_gene_table(
    deg_tables$seurat_layer_adjusted,
    label = "seurat_layer_adjusted"
  )

  ## per-context eQTL overlap remains tied to current Seurat tensorQTL strata.
  deg_by_seurat_context <- lapply(names(SEURAT_CONTEXT_TO_DATASET_ID), function(context_name) {
    col_name <- paste0("n_ttest_sig_", context_name)
    lr_df <- deg_tables$seurat_layer_restricted
    if (!(col_name %in% names(lr_df))) {
      stop("Missing required context localization column in seurat_layer_restricted: ", col_name)
    }

    localized_lr <- lr_df[!is.na(lr_df[[col_name]]) & lr_df[[col_name]] > 0, , drop = FALSE]
    collapse_gene_table(
      rbind(seurat_layer_adjusted_genes, localized_lr[, c("gene_id", "gene_name"), drop = FALSE]),
      label = paste0("deg_by_seurat_context:", context_name)
    )
  })
  names(deg_by_seurat_context) <- names(SEURAT_CONTEXT_TO_DATASET_ID)

  deg_by_dataset_id <- dataset_id_map_from_context_sets(deg_by_seurat_context)
  deg_views <- build_deg_views(
    deg_global = deg_global,
    deg_tables = deg_tables,
    contexts = names(SEURAT_CONTEXT_TO_DATASET_ID)
  )
  disorder_related <- if (isTRUE(include_disorder_degs)) {
    build_disorder_related_degs(
      deg_global = deg_global,
      deg_tables = deg_tables
    )
  } else {
    list(global = list(), validation = list())
  }
  author_union_genes <- read_author_union_list(author_union_file)

  ## author text list validates the reconstructed union; overlaps use tables above.
  author_union_validation <- validate_author_union(deg_global, author_union_genes)

  deg_validation <- list(
    source_counts = vapply(deg_tables, nrow, integer(1)),
    global_gene_id_count = nrow(deg_global),
    global_gene_name_count = length(unique(deg_global$gene_name)),
    seurat_lr_localized_gene_name_counts = vapply(
      names(SEURAT_CONTEXT_TO_DATASET_ID),
      function(context_name) {
        col_name <- paste0("n_ttest_sig_", context_name)
        sum(deg_tables$seurat_layer_restricted[[col_name]] > 0, na.rm = TRUE)
      },
      integer(1)
    ),
    seurat_context_gene_name_counts = vapply(
      deg_by_seurat_context,
      function(df) length(unique(df$gene_name)),
      integer(1)
    ),
    deg_view_gene_counts = deg_views$counts,
    disorder_related = disorder_related$validation,
    author_union = author_union_validation
  )

  files <- c(
    deg_file_paths,
    author_union = author_union_file
  )

  if (isTRUE(verbose)) {
    cat("Loaded author-recommended DEG summaries.\n")
    cat("Repo root:", repo_root, "\n")
    for (nm in names(deg_tables)) {
      cat(sprintf("  %s: %d rows\n", nm, nrow(deg_tables[[nm]])))
    }
    cat("Global union (gene_id):", deg_validation$global_gene_id_count, "\n")
    cat("Global union (gene_name):", deg_validation$global_gene_name_count, "\n")
    if (length(disorder_related$global) > 0) {
      cat("Disorder-related DEG sets (gene_id):\n")
      for (dis in names(disorder_related$global)) {
        cat(sprintf("  %s: %d\n", dis, nrow(disorder_related$global[[dis]])))
      }
    }
    cat("Author union match:", deg_validation$author_union$matches, "\n")
    cat("Per-context Seurat overlap set sizes (gene_name):\n")
    for (ctx in names(deg_validation$seurat_context_gene_name_counts)) {
      ds_id <- unname(SEURAT_CONTEXT_TO_DATASET_ID[[ctx]])
      cat(sprintf("  %s (%s): %d\n", ctx, ds_id, deg_validation$seurat_context_gene_name_counts[[ctx]]))
    }
  }

  list(
    tables = deg_tables,
    global = deg_global,
    by_seurat_context = deg_by_seurat_context,
    by_dataset_id = deg_by_dataset_id,
    disorder_related = disorder_related,
    views = deg_views,
    validation = deg_validation,
    files = files
  )
}

load_custom_cluster_DEGs <- function(repo_root = NULL, verbose = TRUE,
                                     host = JHPCE_HOST, remote_root = JHPCE_REPO_ROOT) {
  repo_root <- resolve_repo_root(repo_root)

  custom_file_paths <- lapply(CUSTOM_CLUSTER_DEG_FILE_SPECS, function(rel_path) {
    ensure_local_file(
      repo_root = repo_root,
      rel_path = rel_path,
      host = host,
      remote_root = remote_root
    )
  })
  custom_tables <- Map(
    f = function(path, label) load_deg_summary(path = path, label = label),
    path = custom_file_paths,
    label = names(custom_file_paths)
  )
  custom_tables <- custom_tables[names(CUSTOM_CLUSTER_DEG_FILE_SPECS)]

  custom_global <- collapse_gene_table(
    do.call(
      rbind,
      lapply(unname(custom_tables), function(df) df[, c("gene_id", "gene_name"), drop = FALSE])
    ),
    label = "custom-cluster DEG summaries"
  )

  custom_layer_adjusted_genes <- collapse_gene_table(
    custom_tables$custom_layer_adjusted,
    label = "custom_layer_adjusted"
  )

  custom_lr_df <- custom_tables$custom_layer_restricted
  custom_by_context <- lapply(names(CUSTOM_CONTEXT_TO_DATASET_ID), function(context_name) {
    col_name <- paste0("n_ttest_sig_", context_name)
    if (!(col_name %in% names(custom_lr_df))) {
      stop("Missing required custom-cluster DEG column: ", col_name)
    }

    localized_lr <- custom_lr_df[!is.na(custom_lr_df[[col_name]]) & custom_lr_df[[col_name]] > 0, , drop = FALSE]
    collapse_gene_table(
      rbind(custom_layer_adjusted_genes, localized_lr[, c("gene_id", "gene_name"), drop = FALSE]),
      label = paste0("deg_by_custom_context:", context_name)
    )
  })
  names(custom_by_context) <- names(CUSTOM_CONTEXT_TO_DATASET_ID)

  custom_by_dataset_id <- dataset_id_map_from_context_sets(
    custom_by_context,
    context_to_dataset_id = CUSTOM_CONTEXT_TO_DATASET_ID
  )
  custom_views <- build_custom_cluster_deg_views(
    custom_by_context = custom_by_context,
    splits = c("all")
  )

  deg_validation <- list(
    source_counts = vapply(custom_tables, nrow, integer(1)),
    global_gene_id_count = nrow(custom_global),
    global_gene_name_count = length(unique(custom_global$gene_name)),
    custom_lr_localized_gene_name_counts = vapply(
      names(CUSTOM_CONTEXT_TO_DATASET_ID),
      function(context_name) {
        col_name <- paste0("n_ttest_sig_", context_name)
        sum(custom_lr_df[[col_name]] > 0, na.rm = TRUE)
      },
      integer(1)
    ),
    custom_context_gene_name_counts = vapply(
      custom_by_context,
      function(df) length(unique(df$gene_name)),
      integer(1)
    ),
    deg_view_gene_counts = custom_views$counts
  )

  files <- custom_file_paths

  if (isTRUE(verbose)) {
    cat("Loaded custom-cluster DEG summaries.\n")
    cat("Repo root:", repo_root, "\n")
    for (nm in names(custom_tables)) {
      cat(sprintf("  %s: %d rows\n", nm, nrow(custom_tables[[nm]])))
    }
    cat("Custom union (gene_id):", deg_validation$global_gene_id_count, "\n")
    cat("Per-context custom overlap set sizes (gene_name):\n")
    for (ctx in names(deg_validation$custom_context_gene_name_counts)) {
      ds_id <- unname(CUSTOM_CONTEXT_TO_DATASET_ID[[ctx]])
      cat(sprintf("  %s (%s): %d\n", ctx, ds_id, deg_validation$custom_context_gene_name_counts[[ctx]]))
    }
  }

  list(
    tables = custom_tables,
    global = custom_global,
    by_custom_context = custom_by_context,
    by_dataset_id = custom_by_dataset_id,
    views = custom_views,
    validation = deg_validation,
    files = files
  )
}

load_DEGs <- function(repo_root = NULL, mode = c("standard", "custom_cluster"), verbose = TRUE,
                      include_disorder_degs = TRUE,
                      host = JHPCE_HOST, remote_root = JHPCE_REPO_ROOT) {
  mode <- match.arg(mode)
  if (identical(mode, "custom_cluster")) {
    return(load_custom_cluster_DEGs(
      repo_root = repo_root,
      verbose = verbose,
      host = host,
      remote_root = remote_root
    ))
  }

  load_standard_DEGs(
    repo_root = repo_root,
    verbose = verbose,
    include_disorder_degs = include_disorder_degs,
    host = host,
    remote_root = remote_root
  )
}

EQTL_REQUIRED_MAP_CIS_COLS <- c("phenotype_id", "variant_id", "qval")
EQTL_REQUIRED_INDEP_COLS <- c("phenotype_id", "variant_id", "rank", "pval_perm")

require_data_table <- function() {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("The data.table package is required")
  }
  invisible(TRUE)
}

path_relative_to <- function(path, root) {
  path_norm <- normalizePath(path, mustWork = FALSE)
  root_norm <- normalizePath(root, mustWork = TRUE)
  prefix <- paste0(root_norm, .Platform$file.sep)
  if (startsWith(path_norm, prefix)) return(sub(paste0("^", prefix), "", path_norm))
  path_norm
}

add_filename_suffix <- function(filename, suffix = "") {
  if (is.null(suffix) || !nzchar(suffix)) return(filename)
  sub("(\\.[^.]+)$", paste0(suffix, "\\1"), filename)
}

load_eqtl_manifest <- function(tqtl_in_dir, context_order = names(SEURAT_CONTEXT_TO_DATASET_ID),
                               split_order = DEG_SPLITS) {
  require_data_table()
  manifest_file <- file.path(tqtl_in_dir, "prep_manifest.csv")
  if (!file.exists(manifest_file)) stop("Missing tensorQTL manifest: ", manifest_file)

  manifest <- data.table::fread(manifest_file)
  req <- c("dataset_id", "seurat_label", "split", "status")
  missing_cols <- setdiff(req, names(manifest))
  if (length(missing_cols) > 0) {
    stop("Missing manifest columns in ", manifest_file, ": ", paste(missing_cols, collapse = ", "))
  }

  if (!"covariate_model" %in% names(manifest)) manifest[, covariate_model := "seurat"]
  manifest <- manifest[
    status == "prepared",
    .(dataset_id, context = seurat_label, split, covariate_model)
  ]
  manifest[, split := data.table::fifelse(
    split == "male",
    "male",
    data.table::fifelse(split == "female", "female", "all")
  )]

  valid_splits <- c("all", "male", "female")
  unexpected_split <- setdiff(unique(manifest$split), valid_splits)
  if (length(unexpected_split) > 0) {
    stop("Unexpected manifest split(s): ", paste(unexpected_split, collapse = ", "))
  }
  manifest <- manifest[split %in% split_order]
  if (nrow(manifest) == 0) {
    stop("No prepared manifest rows match split_order: ", paste(split_order, collapse = ", "))
  }

  if (anyDuplicated(manifest$dataset_id)) {
    stop("Duplicate dataset_id values in manifest: ", manifest_file)
  }
  unexpected_context <- setdiff(unique(manifest$context), context_order)
  if (length(unexpected_context) > 0) {
    stop("Unexpected manifest context(s): ", paste(unexpected_context, collapse = ", "))
  }

  manifest[order(match(split, split_order), match(context, context_order))]
}

read_eqtl_table <- function(dataset_id, tqtl_out_dir, suffix, required_cols,
                            missing_status, repo_root = NULL) {
  require_data_table()
  abs_path <- file.path(tqtl_out_dir, paste0(dataset_id, ".gene.", suffix))
  rel_path <- if (is.null(repo_root)) abs_path else path_relative_to(abs_path, repo_root)

  if (!file.exists(abs_path)) {
    return(list(status = missing_status, rel_path = NA_character_, dt = NULL))
  }
  if (file.info(abs_path)$size == 0) {
    stop("Empty tensorQTL result file: ", rel_path)
  }

  dt <- data.table::fread(abs_path)
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", rel_path, ": ", paste(missing_cols, collapse = ", "))
  }

  list(status = "ok", rel_path = rel_path, dt = dt)
}

read_map_cis <- function(dataset_id, tqtl_out_dir, repo_root = NULL,
                         required_cols = EQTL_REQUIRED_MAP_CIS_COLS) {
  read_eqtl_table(
    dataset_id = dataset_id,
    tqtl_out_dir = tqtl_out_dir,
    suffix = "map_cis.tab.gz",
    required_cols = required_cols,
    missing_status = "missing_map_cis",
    repo_root = repo_root
  )
}

read_map_independent <- function(dataset_id, tqtl_out_dir, repo_root = NULL,
                                 required_cols = EQTL_REQUIRED_INDEP_COLS) {
  read_eqtl_table(
    dataset_id = dataset_id,
    tqtl_out_dir = tqtl_out_dir,
    suffix = "map_independent.txt.gz",
    required_cols = required_cols,
    missing_status = "missing_map_independent",
    repo_root = repo_root
  )
}

collapse_egenes <- function(dt) {
  require_data_table()
  unique(dt[, .(gene_id = phenotype_id)])
}

collapse_gene_list <- function(x) {
  x <- sort(unique(x[!is.na(x) & nzchar(x)]))
  if (length(x) == 0) "" else paste(x, collapse = ", ")
}

normalize_gwas_sets <- function(gwas_sets) {
  if (is.null(gwas_sets)) stop("gwas_sets is required")
  if (!is.list(gwas_sets)) gwas_sets <- list(SCZD = gwas_sets)
  if (is.null(names(gwas_sets)) || any(!nzchar(names(gwas_sets)))) {
    stop("gwas_sets must be named by disorder")
  }

  disorders <- vapply(names(gwas_sets), gwas_check_disorder, character(1))
  names(gwas_sets) <- disorders
  lapply(gwas_sets, function(x) unique(as.character(x[!is.na(x) & nzchar(x)])))
}

getGWASovl <- function(x, dis, gene_col = "gene_name", variant_col = "variant_id",
                       allow_gene_match = TRUE, use_prio = FALSE,
                       gwas_sets = NULL, repo_root = NULL) {
  ## mixed GWASg overlap: variant match OR exact GWAS gene-list match.
  ## variant match means variant_id is in the harmonized significant GWAS set.
  ## gene match means gene_name is in the selected broad GWAS gene list.
  ## set allow_gene_match = FALSE for variant-only GWAS compatibility.
  require_data_table()
  dis <- gwas_check_disorder(dis)

  if (is.data.frame(x)) {
    if (!gene_col %in% names(x)) stop("Missing gene column: ", gene_col)
    gene_name <- as.character(x[[gene_col]])
    variant_id <- if (variant_col %in% names(x)) as.character(x[[variant_col]]) else rep(NA_character_, length(gene_name))
  } else if (is.atomic(x)) {
    gene_name <- as.character(x)
    variant_id <- rep(NA_character_, length(gene_name))
  } else {
    stop("x must be a data.frame or an atomic vector of gene names")
  }

  dt <- data.table::data.table(gene_name = gene_name, variant_id = variant_id)
  gene_list <- loadGWASGeneList(dis, use_prio = use_prio, repo_root = repo_root)
  gene_set <- unique(gene_list$gene_symbol)
  dt[, gene_match := isTRUE(allow_gene_match) & !is.na(gene_name) & nzchar(gene_name) & gene_name %in% gene_set]

  variant_set <- character()
  if (!is.null(gwas_sets)) {
    if (is.list(gwas_sets)) {
      gwas_sets <- normalize_gwas_sets(gwas_sets)
      variant_set <- gwas_sets[[dis]]
    } else {
      variant_set <- unique(as.character(gwas_sets[!is.na(gwas_sets) & nzchar(gwas_sets)]))
    }
  }
  dt[, variant_match := !is.na(variant_id) & nzchar(variant_id) & variant_id %in% variant_set]
  dt[, match_type := data.table::fifelse(
    gene_match & variant_match,
    "gene_and_variant",
    data.table::fifelse(gene_match, "gene", "variant")
  )]
  dt[variant_match == FALSE | is.na(variant_match), variant_id := NA_character_]

  out <- dt[gene_match | variant_match]
  data.table::setcolorder(out, c("gene_name", "variant_id", "gene_match", "variant_match", "match_type"))
  out[]
}

gwas_flag_cols <- function(disorders) paste0(disorders, "_GWAS")
gwasg_flag_cols <- function(disorders) paste0(disorders, "_GWASg")
gwasg_variant_flag_cols <- function(disorders) paste0(gwasg_flag_cols(disorders), "_variant")
gwasg_gene_match_cols <- function(disorders) paste0(gwasg_flag_cols(disorders), "_gene")
gwas_count_cols <- function(disorders) paste0("n_", gwas_flag_cols(disorders))
gwasg_count_cols <- function(disorders) paste0("n_", gwasg_flag_cols(disorders))
deg_gwas_count_cols <- function(disorders) paste0("n_DEG_", gwas_flag_cols(disorders))
deg_gwasg_count_cols <- function(disorders) paste0("n_DEG_", gwasg_flag_cols(disorders))
gwas_gene_cols <- function(disorders) paste0(gwas_flag_cols(disorders), "_genes")
gwasg_gene_cols <- function(disorders) paste0(gwasg_flag_cols(disorders), "_genes")
deg_gwas_gene_cols <- function(disorders) paste0("DEG_", gwas_flag_cols(disorders), "_genes")
deg_gwasg_gene_cols <- function(disorders) paste0("DEG_", gwasg_flag_cols(disorders), "_genes")

paired_gwas_cols <- function(disorders, left_fun, right_fun) {
  as.vector(rbind(left_fun(disorders), right_fun(disorders)))
}

gwas_flag_pair_cols <- function(disorders) paired_gwas_cols(disorders, gwas_flag_cols, gwasg_flag_cols)
gwasg_source_pair_cols <- function(disorders) paired_gwas_cols(disorders, gwasg_variant_flag_cols, gwasg_gene_match_cols)
gwas_count_pair_cols <- function(disorders) paired_gwas_cols(disorders, gwas_count_cols, gwasg_count_cols)
deg_gwas_count_pair_cols <- function(disorders) paired_gwas_cols(disorders, deg_gwas_count_cols, deg_gwasg_count_cols)
gwas_gene_pair_cols <- function(disorders) paired_gwas_cols(disorders, gwas_gene_cols, gwasg_gene_cols)
deg_gwas_gene_pair_cols <- function(disorders) paired_gwas_cols(disorders, deg_gwas_gene_cols, deg_gwasg_gene_cols)

disorder_deg_sets <- function(disorder_related) {
  if (is.null(disorder_related)) return(list())
  if (!is.null(disorder_related$global)) return(disorder_related$global)
  disorder_related
}

disorder_deg_count_cols <- function(disorders) paste0("n_", disorders, "_DEG")
disorder_deg_gene_cols <- function(disorders) paste0(disorders, "_DEG_genes")
disorder_deg_gwas_count_cols <- function(disorders) paste0("n_", disorders, "_DEG_", disorders, "_GWAS")
disorder_deg_gwasg_count_cols <- function(disorders) paste0("n_", disorders, "_DEG_", disorders, "_GWASg")
disorder_deg_gwas_gene_cols <- function(disorders) paste0(disorders, "_DEG_", disorders, "_GWAS_genes")
disorder_deg_gwasg_gene_cols <- function(disorders) paste0(disorders, "_DEG_", disorders, "_GWASg_genes")
disorder_deg_gwas_count_pair_cols <- function(disorders) paired_gwas_cols(
  disorders,
  disorder_deg_gwas_count_cols,
  disorder_deg_gwasg_count_cols
)
disorder_deg_gwas_gene_pair_cols <- function(disorders) paired_gwas_cols(
  disorders,
  disorder_deg_gwas_gene_cols,
  disorder_deg_gwasg_gene_cols
)

order_summary_cols <- function(dt, gwas_disorders, include_signal_count = FALSE) {
  front_cols <- c(
    "split", "context",
    if (include_signal_count) "n_independent_signals",
    "n_eGenes",
    gwas_count_pair_cols(gwas_disorders),
    "n_DEG",
    deg_gwas_count_pair_cols(gwas_disorders),
    gwas_gene_pair_cols(gwas_disorders),
    "DEG_genes",
    deg_gwas_gene_pair_cols(gwas_disorders)
  )
  data.table::setcolorder(dt, c(intersect(front_cols, names(dt)), setdiff(names(dt), front_cols)))
  dt[]
}

order_deg_view_summary_cols <- function(dt, gwas_disorders) {
  front_cols <- c(
    "dataset_id", "context", "split", "deg_view", "deg_sex",
    "n_DEG_reference", "n_eQTL_records", "n_eGenes",
    gwas_count_pair_cols(gwas_disorders),
    "n_eGene_DEG_overlap",
    deg_gwas_count_pair_cols(gwas_disorders),
    "eGene_DEG_overlap_gene_ids", "eGene_DEG_overlap_gene_names",
    deg_gwas_gene_pair_cols(gwas_disorders)
  )
  data.table::setcolorder(dt, c(intersect(front_cols, names(dt)), setdiff(names(dt), front_cols)))
  dt[]
}

gwas_disorders_from_eqtl <- function(dt) {
  disorders <- sub("_GWAS$", "", grep("^[A-Z0-9]+_GWAS$", names(dt), value = TRUE))
  disorders[disorders %in% names(GWAS_BCF_FILES)]
}

gwasx_flag_cols <- function(disorders) paste0(disorders, "_GWASx")
gwasxg_flag_cols <- function(disorders) paste0(disorders, "_GWASxg")
gwasxg_variant_flag_cols <- function(disorders) paste0(gwasxg_flag_cols(disorders), "_variant")
gwasxg_gene_match_cols <- function(disorders) paste0(gwasxg_flag_cols(disorders), "_gene")
gwasx_p_cols <- function(disorders) paste0(disorders, "_GWASx_p")
gwasx_beta_cols <- function(disorders) paste0(disorders, "_GWASx_beta")
gwasx_beta_se_cols <- function(disorders) paste0(disorders, "_GWASx_beta_se")
gwasx_flag_pair_cols <- function(disorders) paired_gwas_cols(disorders, gwasx_flag_cols, gwasxg_flag_cols)
gwasxg_source_pair_cols <- function(disorders) paired_gwas_cols(disorders, gwasxg_variant_flag_cols, gwasxg_gene_match_cols)

gwasx_disorders_from_eqtl <- function(dt) {
  disorders <- sub("_GWASx$", "", grep("^[A-Z0-9]+_GWASx$", names(dt), value = TRUE))
  disorders[disorders %in% names(GWAS_BCF_FILES)]
}

gwas_threshold_from_cache <- function(dis) {
  dis <- gwas_check_disorder(dis)
  if (dis %in% names(GWAS_EXPLORATORY_P_THRESHOLDS)) {
    return(max(GWAS_STRICT_P_THRESHOLD, GWAS_EXPLORATORY_P_THRESHOLDS[[dis]]))
  }
  GWAS_STRICT_P_THRESHOLD
}

load_matched_gwas_by_disorder <- function(disorders = DEFAULT_GWAS_OVERLAP_DISORDERS,
                                          plink2_prefix, repo_root = NULL,
                                          si_min = GWAS_MATCH_SI_MIN) {
  out <- lapply(disorders, function(dis) {
    dis <- gwas_check_disorder(dis)
    ## query broad enough for strict GWAS and any GWASx reuse.
    gwas <- loadGWAS(dis, gwas_threshold_from_cache(dis), repo_root = repo_root, si_min = si_min)
    matched <- matchGwasGeno(gwas, plink2_prefix = plink2_prefix)
    missing_cols <- setdiff(c("variant_id", "p", "beta", "beta_se"), names(matched))
    if (length(missing_cols) > 0) {
      stop("Missing required ", dis, " matched GWAS columns: ", paste(missing_cols, collapse = ", "))
    }
    matched
  })
  names(out) <- vapply(disorders, gwas_check_disorder, character(1))
  out
}

strict_gwas_sets_from_matched <- function(gwas_matched,
                                          p_threshold = GWAS_STRICT_P_THRESHOLD) {
  lapply(gwas_matched, function(dt) {
    dt <- data.table::as.data.table(dt)
    unique(dt[p <= p_threshold & !is.na(variant_id), as.character(variant_id)])
  })
}

gwasx_matched_from_gwas <- function(gwas_matched,
                                    thresholds = GWAS_EXPLORATORY_P_THRESHOLDS) {
  out <- lapply(names(thresholds), function(dis) {
    if (!dis %in% names(gwas_matched)) stop("Missing matched GWAS for GWASx disorder: ", dis)
    dt <- data.table::copy(data.table::as.data.table(gwas_matched[[dis]]))
    ## use a strict suggestive cutoff for exploratory exact-variant matches.
    dt[p < thresholds[[dis]] & !is.na(variant_id), .(
      variant_id = as.character(variant_id),
      p = as.numeric(p),
      beta = as.numeric(beta),
      beta_se = as.numeric(beta_se)
    )]
  })
  names(out) <- names(thresholds)
  out
}

append_disorder_deg_flags <- function(dt, disorder_related,
                                      disorders = DEFAULT_GWASX_DISORDERS) {
  out <- data.table::copy(data.table::as.data.table(dt))
  sets <- disorder_deg_sets(disorder_related)
  for (dis in intersect(disorders, names(sets))) {
    flag_col <- paste0(dis, "_DEG")
    gene_ids <- unique(data.table::as.data.table(sets[[dis]])$gene_id)
    out[, (flag_col) := as.integer(gene_id %in% gene_ids)]
  }
  out[]
}

append_gwasx_stats <- function(dt, gwasx_stats) {
  out <- data.table::copy(data.table::as.data.table(dt))
  for (dis in names(gwasx_stats)) {
    p_col <- paste0(dis, "_GWASx_p")
    beta_col <- paste0(dis, "_GWASx_beta")
    beta_se_col <- paste0(dis, "_GWASx_beta_se")
    out[, (p_col) := NA_real_]
    out[, (beta_col) := NA_real_]
    out[, (beta_se_col) := NA_real_]

    ## p/beta values are defined only for exact variant GWASx matches.
    stats <- data.table::as.data.table(gwasx_stats[[dis]])[, .(
      variant_id,
      gwas_p = p,
      gwas_beta = beta,
      gwas_beta_se = beta_se
    )]
    out[stats, on = "variant_id", (c(p_col, beta_col, beta_se_col)) := list(
      i.gwas_p,
      i.gwas_beta,
      i.gwas_beta_se
    )]
  }

  order_gwasx_cols(out)
}

apply_gwasx_annotations <- function(dt, gwasx_stats, repo_root = NULL) {
  out <- data.table::copy(data.table::as.data.table(dt))
  gwasx_stats <- lapply(gwasx_stats, data.table::as.data.table)
  gwasx_sets <- normalize_gwas_sets(lapply(gwasx_stats, function(x) unique(x$variant_id)))

  for (dis in names(gwasx_sets)) {
    flag_col <- paste0(dis, "_GWASx")
    gwasxg_col <- paste0(dis, "_GWASxg")
    gwasxg_variant_col <- paste0(gwasxg_col, "_variant")
    gwasxg_gene_col <- paste0(gwasxg_col, "_gene")

    out[, (flag_col) := as.integer(variant_id %in% gwasx_sets[[dis]])]
    ovl <- getGWASovl(
      out[, .(gene_name, variant_id)],
      dis = dis,
      gwas_sets = gwasx_sets,
      repo_root = repo_root
    )
    ## GWASxg mirrors GWASg: exact variant match OR broad gene-list match.
    out[, (gwasxg_variant_col) := get(flag_col)]
    out[, (gwasxg_gene_col) := as.integer(!is.na(gene_name) & gene_name %in% ovl[gene_match == TRUE, unique(gene_name)])]
    out[, (gwasxg_col) := as.integer(get(gwasxg_variant_col) == 1L | get(gwasxg_gene_col) == 1L)]
  }

  append_gwasx_stats(out, gwasx_stats)
}

order_gwasx_cols <- function(dt) {
  strict_disorders <- gwas_disorders_from_eqtl(dt)
  gwasx_disorders <- gwasx_disorders_from_eqtl(dt)
  front_cols <- c(
    "source", "result_source", "pair_provenance", "cis_supported", "indep_supported",
    "dataset_id", "context", "split", "gene_id", "gene_name", "DEG",
    paste0(gwasx_disorders, "_DEG"),
    gwas_flag_pair_cols(strict_disorders),
    gwasg_source_pair_cols(strict_disorders),
    gwasx_flag_pair_cols(gwasx_disorders),
    gwasxg_source_pair_cols(gwasx_disorders),
    as.vector(rbind(
      gwasx_p_cols(gwasx_disorders),
      gwasx_beta_cols(gwasx_disorders),
      gwasx_beta_se_cols(gwasx_disorders)
    )),
    "variant_id"
  )
  data.table::setcolorder(dt, c(intersect(front_cols, names(dt)), setdiff(names(dt), front_cols)))
  dt[]
}

summarize_gwasx_significant_pairs <- function(dt, context_order, split_order,
                                              disorder_related = NULL,
                                              gwasx_disorders = gwasx_disorders_from_eqtl(dt)) {
  tmp <- data.table::copy(data.table::as.data.table(dt))
  for (dis in gwasx_disorders) {
    canonical <- c(
      paste0(dis, "_GWAS"),
      paste0(dis, "_GWASg"),
      paste0(dis, "_GWASg_variant"),
      paste0(dis, "_GWASg_gene")
    )
    gwasx <- c(
      paste0(dis, "_GWASx"),
      paste0(dis, "_GWASxg"),
      paste0(dis, "_GWASxg_variant"),
      paste0(dis, "_GWASxg_gene")
    )
    tmp[, (intersect(canonical, names(tmp))) := NULL]
    data.table::setnames(tmp, intersect(gwasx, names(tmp)), canonical[match(intersect(gwasx, names(tmp)), gwasx)])
  }

  out <- summarize_significant_pairs(
    tmp,
    context_order = context_order,
    split_order = split_order,
    disorder_related = disorder_related,
    gwas_disorders = gwasx_disorders
  )
  for (dis in gwasx_disorders) {
    tmp_token <- paste0("__", dis, "_XG_TMP__")
    new_names <- gsub(paste0(dis, "_GWASg"), tmp_token, names(out), fixed = TRUE)
    new_names <- gsub(paste0(dis, "_GWAS"), paste0(dis, "_GWASx"), new_names, fixed = TRUE)
    new_names <- gsub(tmp_token, paste0(dis, "_GWASxg"), new_names, fixed = TRUE)
    data.table::setnames(out, new_names)
  }
  out[]
}

add_gwasx_to_eqtl_result <- function(result, gwasx_stats, degs, repo_root,
                                     context_order, split_order) {
  for (nm in c("map_cis_significant", "map_independent_significant",
               "map_significant_unified", "map_significant_pairs")) {
    if (!is.null(result[[nm]])) {
      result[[nm]] <- append_disorder_deg_flags(
        apply_gwasx_annotations(result[[nm]], gwasx_stats, repo_root = repo_root),
        degs$disorder_related,
        disorders = names(gwasx_stats)
      )
    }
  }
  result$map_significant_summary_GWASx <- summarize_gwasx_significant_pairs(
    result$map_significant_pairs,
    context_order = context_order,
    split_order = split_order,
    disorder_related = degs$disorder_related,
    gwasx_disorders = names(gwasx_stats)
  )
  result
}

tag_eqtl_results <- function(dt, dataset_id, context, split, deg_dt,
                             g2sym, gwas_sets, repo_root = NULL) {
  require_data_table()
  gwas_sets <- normalize_gwas_sets(gwas_sets)
  out <- data.table::copy(dt)
  if ("V1" %in% names(out) && identical(out$V1, seq_len(nrow(out)) - 1L)) {
    out[, V1 := NULL]
  }
  deg_gene_ids <- unique(data.table::as.data.table(deg_dt)$gene_id)

  ## add project annotations while preserving tensorQTL-native columns.
  out[, `:=`(
    dataset_id = dataset_id,
    context = context,
    split = split,
    gene_id = phenotype_id,
    gene_name = g2sym[phenotype_id],
    DEG = as.integer(phenotype_id %in% deg_gene_ids)
  )]
  for (dis in names(gwas_sets)) {
    flag_col <- paste0(dis, "_GWAS")
    gwasg_col <- paste0(dis, "_GWASg")
    gwasg_variant_col <- paste0(gwasg_col, "_variant")
    gwasg_gene_col <- paste0(gwasg_col, "_gene")

    out[, (flag_col) := as.integer(variant_id %in% gwas_sets[[dis]])]
    ovl <- getGWASovl(
      out[, .(gene_name, variant_id)],
      dis = dis,
      gwas_sets = gwas_sets,
      repo_root = repo_root
    )
    out[, (gwasg_variant_col) := get(flag_col)]
    out[, (gwasg_gene_col) := as.integer(!is.na(gene_name) & gene_name %in% ovl[gene_match == TRUE, unique(gene_name)])]
    out[, (gwasg_col) := as.integer(get(gwasg_variant_col) == 1 | get(gwasg_gene_col) == 1)]
  }

  front_cols <- c(
    "dataset_id", "context", "split", "gene_id", "gene_name", "DEG",
    gwas_flag_pair_cols(names(gwas_sets)),
    gwasg_source_pair_cols(names(gwas_sets))
  )
  data.table::setcolorder(out, c(front_cols, setdiff(names(out), front_cols)))
  out[]
}

order_annotated_eqtl_cols <- function(dt) {
  gwas_disorders <- gwas_disorders_from_eqtl(dt)
  front_cols <- c(
    "dataset_id", "context", "split", "gene_id", "gene_name", "DEG",
    gwas_flag_pair_cols(gwas_disorders),
    gwasg_source_pair_cols(gwas_disorders)
  )
  data.table::setcolorder(dt, c(front_cols, setdiff(names(dt), front_cols)))
  dt[]
}

unify_significant_eqtls <- function(map_cis_significant, map_independent_significant) {
  require_data_table()
  pair_key <- c("dataset_id", "context", "split", "gene_id", "gene_name", "variant_id")
  cis <- data.table::copy(data.table::as.data.table(map_cis_significant))
  independent <- data.table::copy(data.table::as.data.table(map_independent_significant))

  missing_cis <- setdiff(pair_key, names(cis))
  if (length(missing_cis) > 0) {
    stop("Missing cis key columns: ", paste(missing_cis, collapse = ", "))
  }
  missing_independent <- setdiff(pair_key, names(independent))
  if (length(missing_independent) > 0) {
    stop("Missing independent key columns: ", paste(missing_independent, collapse = ", "))
  }

  ## keep source rows separate because shared pairs can have different statistics.
  cis[, result_source := "cis"]
  independent[, result_source := "independent"]
  out <- data.table::rbindlist(list(cis, independent), use.names = TRUE, fill = TRUE)
  if (nrow(out) == 0) {
    out[, pair_provenance := character()]
    return(out[])
  }

  pair_sources <- out[, .(
    has_cis = any(result_source == "cis"),
    has_independent = any(result_source == "independent")
  ), by = pair_key]
  pair_sources[, pair_provenance := data.table::fcase(
    has_cis & has_independent, "cis_and_independent",
    has_cis, "cis_only",
    has_independent, "independent_only"
  )]

  out[pair_sources[, c(pair_key, "pair_provenance"), with = FALSE],
      pair_provenance := i.pair_provenance,
      on = pair_key]
  front_cols <- c("result_source", "pair_provenance", pair_key)
  data.table::setcolorder(out, c(front_cols, setdiff(names(out), front_cols)))
  out[]
}

collapse_significant_eqtl_pairs <- function(map_significant_unified) {
  require_data_table()
  pair_key <- c("dataset_id", "context", "split", "gene_id", "gene_name", "variant_id")
  dt <- data.table::copy(data.table::as.data.table(map_significant_unified))

  missing_cols <- setdiff(c("result_source", pair_key), names(dt))
  if (length(missing_cols) > 0) {
    stop("Missing unified eQTL columns: ", paste(missing_cols, collapse = ", "))
  }
  bad_sources <- setdiff(unique(dt$result_source), c("cis", "independent"))
  if (length(bad_sources) > 0) {
    stop("Unexpected result_source values: ", paste(bad_sources, collapse = ", "))
  }
  source_key <- c("result_source", pair_key)
  if (anyDuplicated(dt[, source_key, with = FALSE]) > 0) {
    stop("Duplicate source-specific eQTL pair rows found")
  }

  pair_sources <- dt[, .(
    cis_supported = any(result_source == "cis"),
    indep_supported = any(result_source == "independent")
  ), by = pair_key]
  pair_sources[, pair_provenance := data.table::fcase(
    cis_supported & indep_supported, "cis_and_independent",
    cis_supported, "cis_only",
    indep_supported, "independent_only"
  )]

  ## keep canonical cis rows for shared pairs; keep independent-only rows.
  dt[, source_priority := data.table::fifelse(result_source == "cis", 1L, 2L)]
  data.table::setorderv(dt, c(pair_key, "source_priority"))
  out <- dt[, .SD[1], by = pair_key]
  out[, source := data.table::fifelse(result_source == "independent", "indep", "cis")]
  out[pair_sources, `:=`(
    pair_provenance = i.pair_provenance,
    cis_supported = i.cis_supported,
    indep_supported = i.indep_supported
  ), on = pair_key]
  drop_cols <- intersect(c("result_source", "source_priority"), names(out))
  out[, (drop_cols) := NULL]

  front_cols <- c("source", "pair_provenance", "cis_supported", "indep_supported", pair_key)
  data.table::setcolorder(out, c(front_cols, setdiff(names(out), front_cols)))
  out[]
}

complete_summary_grid <- function(context_order = names(SEURAT_CONTEXT_TO_DATASET_ID),
                                  split_order = DEG_SPLITS) {
  require_data_table()
  data.table::CJ(split = split_order, context = context_order, unique = TRUE)
}

fill_summary_missing <- function(out, count_cols, text_cols) {
  for (col in count_cols) data.table::set(out, which(is.na(out[[col]])), col, 0L)
  for (col in text_cols) data.table::set(out, which(is.na(out[[col]])), col, "")
  out
}

summary_template <- function(context_order = names(SEURAT_CONTEXT_TO_DATASET_ID),
                             split_order = DEG_SPLITS,
                             include_signal_count = FALSE,
                             gwas_disorders = "SCZD") {
  require_data_table()
  out <- complete_summary_grid(context_order = context_order, split_order = split_order)
  if (isTRUE(include_signal_count)) out[, n_independent_signals := 0L]
  out[, n_eGenes := 0L]
  for (col in c(
    gwas_count_cols(gwas_disorders), gwasg_count_cols(gwas_disorders),
    deg_gwas_count_cols(gwas_disorders), deg_gwasg_count_cols(gwas_disorders)
  )) {
    out[, (col) := 0L]
  }
  out[, n_DEG := 0L]
  for (col in c(
    gwas_gene_cols(gwas_disorders), gwasg_gene_cols(gwas_disorders),
    deg_gwas_gene_cols(gwas_disorders), deg_gwasg_gene_cols(gwas_disorders)
  )) {
    out[, (col) := ""]
  }
  out[, DEG_genes := ""]
  order_summary_cols(out, gwas_disorders = gwas_disorders, include_signal_count = include_signal_count)
}

summarize_tagged_eqtls <- function(dt, context_order = names(SEURAT_CONTEXT_TO_DATASET_ID),
                                   split_order = DEG_SPLITS,
                                   signal_count = FALSE,
                                   gwas_disorders = gwas_disorders_from_eqtl(dt)) {
  require_data_table()
  include_signal_count <- isTRUE(signal_count)
  if (length(gwas_disorders) == 0) gwas_disorders <- "SCZD"
  if (nrow(dt) == 0) {
    return(summary_template(
      context_order = context_order,
      split_order = split_order,
      include_signal_count = include_signal_count,
      gwas_disorders = gwas_disorders
    ))
  }

  if (include_signal_count) {
    obs <- dt[, .(
      n_independent_signals = .N,
      n_eGenes = data.table::uniqueN(gene_id),
      n_DEG = data.table::uniqueN(gene_id[DEG == 1]),
      DEG_genes = collapse_gene_list(gene_name[DEG == 1])
    ), by = .(split, context)]
    for (dis in gwas_disorders) {
      flag_col <- paste0(dis, "_GWAS")
      gwasg_col <- paste0(dis, "_GWASg")
      obs[, paste0("n_", flag_col) := dt[obs, data.table::uniqueN(gene_id[get(flag_col) == 1]),
        by = .EACHI, on = .(split, context)]$V1]
      obs[, paste0("n_", gwasg_col) := dt[obs, data.table::uniqueN(gene_id[get(gwasg_col) == 1]),
        by = .EACHI, on = .(split, context)]$V1]
      obs[, paste0("n_DEG_", flag_col) := dt[obs, data.table::uniqueN(gene_id[DEG == 1 & get(flag_col) == 1]),
        by = .EACHI, on = .(split, context)]$V1]
      obs[, paste0("n_DEG_", gwasg_col) := dt[obs, data.table::uniqueN(gene_id[DEG == 1 & get(gwasg_col) == 1]),
        by = .EACHI, on = .(split, context)]$V1]
      obs[, paste0(flag_col, "_genes") := dt[obs, collapse_gene_list(gene_name[get(flag_col) == 1]),
        by = .EACHI, on = .(split, context)]$V1]
      obs[, paste0(gwasg_col, "_genes") := dt[obs, collapse_gene_list(gene_name[get(gwasg_col) == 1]),
        by = .EACHI, on = .(split, context)]$V1]
      obs[, paste0("DEG_", flag_col, "_genes") := dt[obs, collapse_gene_list(gene_name[DEG == 1 & get(flag_col) == 1]),
        by = .EACHI, on = .(split, context)]$V1]
      obs[, paste0("DEG_", gwasg_col, "_genes") := dt[obs, collapse_gene_list(gene_name[DEG == 1 & get(gwasg_col) == 1]),
        by = .EACHI, on = .(split, context)]$V1]
    }
  } else {
    obs <- dt[, .(
      n_eGenes = data.table::uniqueN(gene_id),
      n_DEG = data.table::uniqueN(gene_id[DEG == 1]),
      DEG_genes = collapse_gene_list(gene_name[DEG == 1])
    ), by = .(split, context)]
    for (dis in gwas_disorders) {
      flag_col <- paste0(dis, "_GWAS")
      gwasg_col <- paste0(dis, "_GWASg")
      obs[, paste0("n_", flag_col) := dt[obs, data.table::uniqueN(gene_id[get(flag_col) == 1]),
        by = .EACHI, on = .(split, context)]$V1]
      obs[, paste0("n_", gwasg_col) := dt[obs, data.table::uniqueN(gene_id[get(gwasg_col) == 1]),
        by = .EACHI, on = .(split, context)]$V1]
      obs[, paste0("n_DEG_", flag_col) := dt[obs, data.table::uniqueN(gene_id[DEG == 1 & get(flag_col) == 1]),
        by = .EACHI, on = .(split, context)]$V1]
      obs[, paste0("n_DEG_", gwasg_col) := dt[obs, data.table::uniqueN(gene_id[DEG == 1 & get(gwasg_col) == 1]),
        by = .EACHI, on = .(split, context)]$V1]
      obs[, paste0(flag_col, "_genes") := dt[obs, collapse_gene_list(gene_name[get(flag_col) == 1]),
        by = .EACHI, on = .(split, context)]$V1]
      obs[, paste0(gwasg_col, "_genes") := dt[obs, collapse_gene_list(gene_name[get(gwasg_col) == 1]),
        by = .EACHI, on = .(split, context)]$V1]
      obs[, paste0("DEG_", flag_col, "_genes") := dt[obs, collapse_gene_list(gene_name[DEG == 1 & get(flag_col) == 1]),
        by = .EACHI, on = .(split, context)]$V1]
      obs[, paste0("DEG_", gwasg_col, "_genes") := dt[obs, collapse_gene_list(gene_name[DEG == 1 & get(gwasg_col) == 1]),
        by = .EACHI, on = .(split, context)]$V1]
    }
  }

  out <- merge(
    summary_template(
      context_order = context_order,
      split_order = split_order,
      include_signal_count = include_signal_count,
      gwas_disorders = gwas_disorders
    )[, .(split, context)],
    obs,
    by = c("split", "context"),
    all.x = TRUE
  )
  count_cols <- c(
    "n_eGenes", gwas_count_cols(gwas_disorders), gwasg_count_cols(gwas_disorders),
    "n_DEG", deg_gwas_count_cols(gwas_disorders), deg_gwasg_count_cols(gwas_disorders)
  )
  if (include_signal_count) count_cols <- c("n_independent_signals", count_cols)
  fill_summary_missing(
    out,
    count_cols = count_cols,
    text_cols = c(
      gwas_gene_cols(gwas_disorders), gwasg_gene_cols(gwas_disorders),
      "DEG_genes", deg_gwas_gene_cols(gwas_disorders), deg_gwasg_gene_cols(gwas_disorders)
    )
  )
  out <- order_summary_cols(
    out[order(match(split, split_order), match(context, context_order))],
    gwas_disorders = gwas_disorders,
    include_signal_count = include_signal_count
  )
  out
}

summarize_independent_signals <- function(dt, context_order = names(SEURAT_CONTEXT_TO_DATASET_ID),
                                          split_order = DEG_SPLITS,
                                          gwas_disorders = gwas_disorders_from_eqtl(dt)) {
  summarize_tagged_eqtls(
    dt,
    context_order = context_order,
    split_order = split_order,
    signal_count = TRUE,
    gwas_disorders = gwas_disorders
  )
}

summarize_disorder_deg_overlaps <- function(dt, disorder_related,
                                            context_order = names(SEURAT_CONTEXT_TO_DATASET_ID),
                                            split_order = DEG_SPLITS,
                                            gwas_disorders = gwas_disorders_from_eqtl(dt)) {
  require_data_table()
  sets <- disorder_deg_sets(disorder_related)
  if (length(sets) == 0) {
    return(complete_summary_grid(context_order = context_order, split_order = split_order))
  }

  out <- complete_summary_grid(context_order = context_order, split_order = split_order)
  work <- data.table::copy(data.table::as.data.table(dt))
  disorder_names <- names(sets)
  paired_disorders <- intersect(disorder_names, gwas_disorders)

  for (dis in disorder_names) {
    flag_col <- paste0(dis, "_DEG")
    count_col <- paste0("n_", flag_col)
    genes_col <- paste0(flag_col, "_genes")
    gene_ids <- unique(data.table::as.data.table(sets[[dis]])$gene_id)

    ## keep disorder DEG flags internal to summary construction.
    work[, (flag_col) := as.integer(gene_id %in% gene_ids)]
    obs <- work[, .(
      count = data.table::uniqueN(gene_id[get(flag_col) == 1]),
      genes = collapse_gene_list(gene_name[get(flag_col) == 1])
    ), by = .(split, context)]
    data.table::setnames(obs, c("count", "genes"), c(count_col, genes_col))
    out <- merge(out, obs, by = c("split", "context"), all.x = TRUE)

    if (dis %in% paired_disorders) {
      gwas_col <- paste0(dis, "_GWAS")
      gwasg_col <- paste0(dis, "_GWASg")
      gwas_count_col <- paste0("n_", flag_col, "_", gwas_col)
      gwasg_count_col <- paste0("n_", flag_col, "_", gwasg_col)
      gwas_genes_col <- paste0(flag_col, "_", gwas_col, "_genes")
      gwasg_genes_col <- paste0(flag_col, "_", gwasg_col, "_genes")

      obs_gwas <- work[, .(
        gwas_count = data.table::uniqueN(gene_id[get(flag_col) == 1 & get(gwas_col) == 1]),
        gwasg_count = data.table::uniqueN(gene_id[get(flag_col) == 1 & get(gwasg_col) == 1]),
        gwas_genes = collapse_gene_list(gene_name[get(flag_col) == 1 & get(gwas_col) == 1]),
        gwasg_genes = collapse_gene_list(gene_name[get(flag_col) == 1 & get(gwasg_col) == 1])
      ), by = .(split, context)]
      data.table::setnames(
        obs_gwas,
        c("gwas_count", "gwasg_count", "gwas_genes", "gwasg_genes"),
        c(gwas_count_col, gwasg_count_col, gwas_genes_col, gwasg_genes_col)
      )
      out <- merge(out, obs_gwas, by = c("split", "context"), all.x = TRUE)
    }
  }

  count_cols <- c(
    disorder_deg_count_cols(disorder_names),
    disorder_deg_gwas_count_pair_cols(paired_disorders)
  )
  text_cols <- c(
    disorder_deg_gene_cols(disorder_names),
    disorder_deg_gwas_gene_pair_cols(paired_disorders)
  )
  fill_summary_missing(out, count_cols = count_cols, text_cols = text_cols)
}

summarize_significant_pairs <- function(dt, context_order = names(SEURAT_CONTEXT_TO_DATASET_ID),
                                        split_order = DEG_SPLITS,
                                        disorder_related = NULL,
                                        gwas_disorders = gwas_disorders_from_eqtl(dt)) {
  require_data_table()
  if (!"source" %in% names(dt)) stop("Missing source column in significant pair table")
  req <- c("pair_provenance", "cis_supported", "indep_supported")
  missing_cols <- setdiff(req, names(dt))
  if (length(missing_cols) > 0) {
    stop("Missing significant pair support columns: ", paste(missing_cols, collapse = ", "))
  }
  if (length(gwas_disorders) == 0) gwas_disorders <- "SCZD"

  out <- summarize_tagged_eqtls(
    dt,
    context_order = context_order,
    split_order = split_order,
    signal_count = FALSE,
    gwas_disorders = gwas_disorders
  )
  disorder_summary <- summarize_disorder_deg_overlaps(
    dt = dt,
    disorder_related = disorder_related,
    context_order = context_order,
    split_order = split_order,
    gwas_disorders = gwas_disorders
  )
  out <- merge(out, disorder_summary, by = c("split", "context"), all.x = TRUE)
  if (nrow(dt) > 0) {
    pair_counts <- dt[, .(
      n_significant_pairs = .N,
      n_cis_supported_pairs = sum(cis_supported),
      n_indep_supported_pairs = sum(indep_supported),
      n_shared_pairs = sum(pair_provenance == "cis_and_independent"),
      n_cis_only_pairs = sum(pair_provenance == "cis_only"),
      n_indep_only_pairs = sum(pair_provenance == "independent_only")
    ), by = .(split, context)]
  } else {
    pair_counts <- data.table::data.table(
      split = character(),
      context = character(),
      n_significant_pairs = integer(),
      n_cis_supported_pairs = integer(),
      n_indep_supported_pairs = integer(),
      n_shared_pairs = integer(),
      n_cis_only_pairs = integer(),
      n_indep_only_pairs = integer()
    )
  }
  out <- merge(out, pair_counts, by = c("split", "context"), all.x = TRUE)
  pair_count_cols <- c(
    "n_significant_pairs", "n_cis_supported_pairs", "n_indep_supported_pairs",
    "n_shared_pairs", "n_cis_only_pairs", "n_indep_only_pairs"
  )
  for (col in pair_count_cols) {
    data.table::set(out, which(is.na(out[[col]])), col, 0L)
  }
  disorder_names <- names(disorder_deg_sets(disorder_related))
  paired_disorders <- intersect(disorder_names, gwas_disorders)
  front_cols <- c(
    "split", "context", pair_count_cols,
    "n_eGenes", gwas_count_pair_cols(gwas_disorders),
    "n_DEG", disorder_deg_count_cols(disorder_names),
    deg_gwas_count_pair_cols(gwas_disorders),
    disorder_deg_gwas_count_pair_cols(paired_disorders),
    gwas_gene_pair_cols(gwas_disorders),
    "DEG_genes", disorder_deg_gene_cols(disorder_names),
    deg_gwas_gene_pair_cols(gwas_disorders),
    disorder_deg_gwas_gene_pair_cols(paired_disorders)
  )
  data.table::setcolorder(out, c(front_cols, setdiff(names(out), front_cols)))
  out[order(match(split, split_order), match(context, context_order))]
}

## ---- output column finalization (clean GWAS/DEG overlap naming) ----------
## The internal engine emits legacy overlap columns per disorder DIS:
##   DIS_GWAS / DIS_GWASx          variant in GWAS set at 5e-8 / 1e-5
##   DIS_GWASg_gene (= DIS_GWASxg_gene)  eGene in curated disorder gene list
##   DIS_GWASg / DIS_GWASxg        variant OR gene-list (combined)
##   DIS_GWASg_variant (= DIS_GWAS), DIS_GWASxg_variant (= DIS_GWASx)  duplicates
##   DIS_GWASx_p/_beta/_beta_se    GWAS stats at the matched variant
## These helpers relabel the *output* tables to the cleaned scheme
##   DIS_gwasVar_strict / DIS_gwasVar_exp / DIS_gwasGene
##   DIS_gwas_strict / DIS_gwas_exp / DIS_gwasP / DIS_gwasBeta / DIS_gwasBetaSE
## dropping the redundant *_variant duplicates and the duplicate exp gene column.

canonical_overlap_disorder <- function(dis) {
  tryCatch(gwas_check_disorder(dis), error = function(e) toupper(dis))
}

rename_or_drop_existing_col <- function(dt, old, new) {
  if (!old %in% names(dt)) return(invisible(dt))
  if (identical(old, new)) return(invisible(dt))
  if (new %in% names(dt)) {
    if (!identical(dt[[old]], dt[[new]])) {
      stop("Cannot rename ", old, " to existing column ", new, ": values differ.")
    }
    dt[, (old) := NULL]
    return(invisible(dt))
  }
  data.table::setnames(dt, old, new)
  invisible(dt)
}

eqtl_overlap_disorders <- function(dt) {
  unique(sub("_GWAS.*$", "", grep("^[A-Z0-9]+_GWAS", names(dt), value = TRUE)))
}

order_eqtl_overlap_cols <- function(dt) {
  disorders <- intersect(
    c("MDD", "BD", "SCZD"),
    sub("_gwasVar_strict$", "", grep("_gwasVar_strict$", names(dt), value = TRUE))
  )
  overlap <- unlist(lapply(disorders, function(d) paste0(d, c(
    "_gwasVar_strict", "_gwasVar_exp", "_gwasGene",
    "_gwas_strict", "_gwas_exp", "_gwasP", "_gwasBeta", "_gwasBetaSE"
  ))))
  front <- c(
    intersect(c("source", "result_source", "pair_provenance", "cis_supported",
                "indep_supported", "dataset_id", "context", "split",
                "gene_id", "gene_name", "DEG", "MDD_DEG", "BD_DEG"), names(dt)),
    intersect(overlap, names(dt))
  )
  data.table::setcolorder(dt, c(front, setdiff(names(dt), front)))
  dt[]
}

## relabel per-eQTL-row overlap columns; tolerant of strict-only tables.
rename_eqtl_overlap_cols <- function(dt) {
  require_data_table()
  dt <- data.table::as.data.table(dt)
  for (dis in eqtl_overlap_disorders(dt)) {
    out_dis <- canonical_overlap_disorder(dis)
    drop <- intersect(paste0(dis, c("_GWASg_variant", "_GWASxg_variant", "_GWASxg_gene")), names(dt))
    if (length(drop)) dt[, (drop) := NULL]
    map <- c(
      "_GWAS" = "_gwasVar_strict",
      "_GWASx" = "_gwasVar_exp",
      "_GWASg_gene" = "_gwasGene",
      "_GWASg" = "_gwas_strict",
      "_GWASxg" = "_gwas_exp",
      "_GWASx_p" = "_gwasP",
      "_GWASx_beta" = "_gwasBeta",
      "_GWASx_beta_se" = "_gwasBetaSE"
    )
    old <- paste0(dis, names(map))
    new <- paste0(out_dis, unname(map))
    keep <- old %in% names(dt)
    if (any(keep)) {
      for (i in which(keep)) rename_or_drop_existing_col(dt, old[i], new[i])
    }
  }
  order_eqtl_overlap_cols(dt)
}

## relabel summary overlap columns: drop variant-only counts, keep the combined
## (variant-or-gene-list) counts as gwas_<level>, broad-DEG trifecta as
## trifecta_<DIS>_<level>, and disorder-DEG overlap as <DIS>_DEGxGWAS_<level>.
rename_summary_overlap_cols <- function(dt) {
  require_data_table()
  dt <- data.table::as.data.table(dt)
  relabel_one <- function(dt, dis, combined, variant_only) {
    out_dis <- canonical_overlap_disorder(dis)
    suffix <- if (identical(combined, "GWASg")) "strict" else "exp"
    drop <- intersect(c(
      paste0("n_", dis, "_", variant_only),
      paste0("n_DEG_", dis, "_", variant_only),
      paste0("n_", dis, "_DEG_", dis, "_", variant_only),
      paste0(dis, "_", variant_only, "_genes"),
      paste0("DEG_", dis, "_", variant_only, "_genes"),
      paste0(dis, "_DEG_", dis, "_", variant_only, "_genes")
    ), names(dt))
    if (length(drop)) dt[, (drop) := NULL]
    pairs <- c(
      paste0("n_", dis, "_", combined),                   paste0("n_", out_dis, "_gwas_", suffix),
      paste0(dis, "_", combined, "_genes"),               paste0(out_dis, "_gwas_", suffix, "_genes"),
      paste0("n_DEG_", dis, "_", combined),               paste0("n_trifecta_", out_dis, "_", suffix),
      paste0("DEG_", dis, "_", combined, "_genes"),       paste0("trifecta_", out_dis, "_", suffix, "_genes"),
      paste0("n_", dis, "_DEG_", dis, "_", combined),     paste0("n_", out_dis, "_DEGxGWAS_", suffix),
      paste0(dis, "_DEG_", dis, "_", combined, "_genes"), paste0(out_dis, "_DEGxGWAS_", suffix, "_genes")
    )
    old <- pairs[c(TRUE, FALSE)]
    new <- pairs[c(FALSE, TRUE)]
    keep <- old %in% names(dt)
    if (any(keep)) {
      for (i in which(keep)) rename_or_drop_existing_col(dt, old[i], new[i])
    }
    dt
  }
  for (dis in sub("^n_", "", sub("_GWASg$", "", grep("^n_[A-Z0-9]+_GWASg$", names(dt), value = TRUE)))) {
    dt <- relabel_one(dt, dis, "GWASg", "GWAS")
  }
  for (dis in sub("^n_", "", sub("_GWASxg$", "", grep("^n_[A-Z0-9]+_GWASxg$", names(dt), value = TRUE)))) {
    dt <- relabel_one(dt, dis, "GWASxg", "GWASx")
  }
  dt
}

order_significant_summary_cols <- function(dt) {
  disorders <- intersect(
    c("MDD", "BD", "SCZD"),
    sub("^n_", "", sub("_gwas_strict$", "", grep("^n_[A-Z0-9]+_gwas_strict$", names(dt), value = TRUE)))
  )
  by_dis <- function(tpl) unlist(lapply(disorders, function(d) sprintf(tpl, d, c("strict", "exp"))))
  front <- c(
    "split", "context",
    intersect(c("n_eQTL_signals", "n_signal_member_pairs", "n_independent_signals",
                "n_significant_pairs", "n_cis_supported_pairs",
                "n_indep_supported_pairs", "n_shared_pairs", "n_cis_only_pairs",
                "n_indep_only_pairs", "n_nominal_pairs"), names(dt)),
    "n_eGenes", "n_DEG",
    intersect(c("n_MDD_DEG", "n_BD_DEG"), names(dt)),
    by_dis("n_%s_gwas_%s"), by_dis("n_trifecta_%s_%s"), by_dis("n_%s_DEGxGWAS_%s"),
    "DEG_genes",
    intersect(c("MDD_DEG_genes", "BD_DEG_genes"), names(dt)),
    by_dis("%s_gwas_%s_genes"), by_dis("trifecta_%s_%s_genes"), by_dis("%s_DEGxGWAS_%s_genes")
  )
  front <- intersect(front, names(dt))
  data.table::setcolorder(dt, c(front, setdiff(names(dt), front)))
  dt[]
}

## merge the strict (5e-8) summary with the exploratory (1e-5) summary into a
## single relabelled summary carrying both levels side by side.
build_merged_significant_summary <- function(strict, exp = NULL) {
  require_data_table()
  key <- c("split", "context")
  out <- rename_summary_overlap_cols(data.table::copy(data.table::as.data.table(strict)))
  if (!is.null(exp)) {
    x <- rename_summary_overlap_cols(data.table::copy(data.table::as.data.table(exp)))
    exp_cols <- setdiff(grep("_exp(_genes)?$", names(x), value = TRUE), names(out))
    if (length(exp_cols)) {
      out <- merge(out, x[, c(key, exp_cols), with = FALSE], by = key, all.x = TRUE)
    }
  }
  order_significant_summary_cols(out)
}

summarize_eqtl_deg_views <- function(eqtl_dt, manifest, degs,
                                     context_order = names(SEURAT_CONTEXT_TO_DATASET_ID),
                                     split_order = DEG_SPLITS,
                                     signal_count = FALSE,
                                     gwas_disorders = gwas_disorders_from_eqtl(eqtl_dt)) {
  require_data_table()
  if (length(gwas_disorders) == 0) gwas_disorders <- "SCZD"
  view_dt <- data.table::as.data.table(degs$views$long)
  view_dt <- view_dt[, .(
    deg_view, context, split, deg_sex, gene_id,
    deg_gene_name = gene_name
  )]

  view_counts <- view_dt[, .(
    n_DEG_reference = data.table::uniqueN(gene_id)
  ), by = .(deg_view, context, split, deg_sex)]

  template <- merge(
    manifest[, .(dataset_id, context, split)],
    view_counts,
    by = c("context", "split"),
    all.x = FALSE,
    all.y = TRUE,
    allow.cartesian = TRUE
  )

  totals <- eqtl_dt[, .(
    n_eQTL_records = .N,
    n_eGenes = data.table::uniqueN(phenotype_id)
  ), by = .(dataset_id, context, split)]
  for (dis in gwas_disorders) {
    flag_col <- paste0(dis, "_GWAS")
    gwasg_col <- paste0(dis, "_GWASg")
    out_col <- paste0("n_", flag_col)
    gwasg_out_col <- paste0("n_", gwasg_col)
    totals[, (out_col) := eqtl_dt[totals, data.table::uniqueN(phenotype_id[get(flag_col) == 1]),
      by = .EACHI, on = .(dataset_id, context, split)]$V1]
    totals[, (gwasg_out_col) := eqtl_dt[totals, data.table::uniqueN(phenotype_id[get(gwasg_col) == 1]),
      by = .EACHI, on = .(dataset_id, context, split)]$V1]
  }

  overlap_dt <- merge(
    eqtl_dt,
    view_dt,
    by = c("context", "split", "gene_id"),
    all = FALSE,
    allow.cartesian = TRUE
  )

  overlaps <- overlap_dt[, .(
    n_eGene_DEG_overlap = data.table::uniqueN(gene_id),
    eGene_DEG_overlap_gene_ids = collapse_gene_list(gene_id),
    eGene_DEG_overlap_gene_names = collapse_gene_list(deg_gene_name)
  ), by = .(dataset_id, context, split, deg_view, deg_sex)]
  for (dis in gwas_disorders) {
    flag_col <- paste0(dis, "_GWAS")
    gwasg_col <- paste0(dis, "_GWASg")
    count_col <- paste0("n_DEG_", flag_col)
    gwasg_count_col <- paste0("n_DEG_", gwasg_col)
    genes_col <- paste0("DEG_", flag_col, "_genes")
    gwasg_genes_col <- paste0("DEG_", gwasg_col, "_genes")
    overlaps[, (count_col) := overlap_dt[overlaps, data.table::uniqueN(gene_id[get(flag_col) == 1]),
      by = .EACHI, on = .(dataset_id, context, split, deg_view, deg_sex)]$V1]
    overlaps[, (gwasg_count_col) := overlap_dt[overlaps, data.table::uniqueN(gene_id[get(gwasg_col) == 1]),
      by = .EACHI, on = .(dataset_id, context, split, deg_view, deg_sex)]$V1]
    overlaps[, (genes_col) := overlap_dt[overlaps, collapse_gene_list(deg_gene_name[get(flag_col) == 1]),
      by = .EACHI, on = .(dataset_id, context, split, deg_view, deg_sex)]$V1]
    overlaps[, (gwasg_genes_col) := overlap_dt[overlaps, collapse_gene_list(deg_gene_name[get(gwasg_col) == 1]),
      by = .EACHI, on = .(dataset_id, context, split, deg_view, deg_sex)]$V1]
  }

  out <- merge(template, totals, by = c("dataset_id", "context", "split"), all.x = TRUE)
  out <- merge(
    out,
    overlaps,
    by = c("dataset_id", "context", "split", "deg_view", "deg_sex"),
    all.x = TRUE
  )

  fill_summary_missing(
    out,
    count_cols = c(
      "n_eQTL_records", "n_eGenes",
      gwas_count_cols(gwas_disorders), gwasg_count_cols(gwas_disorders),
      "n_eGene_DEG_overlap",
      deg_gwas_count_cols(gwas_disorders), deg_gwasg_count_cols(gwas_disorders)
    ),
    text_cols = c(
      "eGene_DEG_overlap_gene_ids", "eGene_DEG_overlap_gene_names",
      deg_gwas_gene_cols(gwas_disorders), deg_gwasg_gene_cols(gwas_disorders)
    )
  )

  view_order <- c(
    "custom_primary",
    "broad_interaction",
    "context_localized",
    "sex_specific",
    "context_and_sex_specific"
  )
  out <- out[order(
    match(split, split_order),
    match(context, context_order),
    match(deg_view, view_order),
    deg_sex
  )]
  order_deg_view_summary_cols(out, gwas_disorders = gwas_disorders)
}

summarize_dataset_counts <- function(dataset_id, context, split, file_info) {
  require_data_table()
  if (identical(file_info$status, "missing_map_independent")) {
    return(data.table::data.table(
      dataset_id = dataset_id,
      context = context,
      split = split,
      status = file_info$status,
      map_independent_file = NA_character_,
      n_independent_signals = NA_integer_,
      n_eGenes = NA_integer_,
      max_rank = NA_integer_,
      lead_variant_count = NA_integer_
    ))
  }

  dt <- file_info$dt
  data.table::data.table(
    dataset_id = dataset_id,
    context = context,
    split = split,
    status = file_info$status,
    map_independent_file = file_info$rel_path,
    n_independent_signals = nrow(dt),
    n_eGenes = data.table::uniqueN(dt$phenotype_id),
    max_rank = as.integer(max(dt$rank, na.rm = TRUE)),
    lead_variant_count = data.table::uniqueN(dt$variant_id)
  )
}

summarize_dataset_overlap <- function(dataset_id, context, split, file_info, deg_dt) {
  require_data_table()
  n_context_degs <- nrow(deg_dt)

  if (identical(file_info$status, "missing_map_independent")) {
    return(data.table::data.table(
      dataset_id = dataset_id,
      context = context,
      split = split,
      status = file_info$status,
      map_independent_file = NA_character_,
      n_independent_signals = NA_integer_,
      n_eGenes = NA_integer_,
      n_context_DEGs = n_context_degs,
      n_eGene_DEG_overlap = NA_integer_,
      eGene_DEG_overlap_gene_ids = NA_character_,
      eGene_DEG_overlap_gene_names = NA_character_
    ))
  }

  dt <- file_info$dt
  egenes <- collapse_egenes(dt)
  overlap_dt <- merge(
    egenes,
    data.table::as.data.table(deg_dt),
    by = "gene_id",
    all = FALSE,
    sort = FALSE
  )

  data.table::data.table(
    dataset_id = dataset_id,
    context = context,
    split = split,
    status = file_info$status,
    map_independent_file = file_info$rel_path,
    n_independent_signals = nrow(dt),
    n_eGenes = nrow(egenes),
    n_context_DEGs = n_context_degs,
    n_eGene_DEG_overlap = data.table::uniqueN(overlap_dt$gene_id),
    eGene_DEG_overlap_gene_ids = collapse_gene_list(overlap_dt$gene_id),
    eGene_DEG_overlap_gene_names = collapse_gene_list(overlap_dt$gene_name)
  )
}

manifest_qc_table <- function(tqtl_in_dir, split_order = DEG_SPLITS) {
  require_data_table()
  manifest_file <- file.path(tqtl_in_dir, "prep_manifest.csv")
  if (!file.exists(manifest_file)) stop("Missing tensorQTL manifest: ", manifest_file)
  manifest <- data.table::fread(manifest_file)
  keep_cols <- intersect(
    c(
      "dataset_id", "seurat_label", "cluster_col", "analysis_label", "split",
      "covariate_model", "n_samples", "n_genes_bed", "n_genes_pca",
      "n_expr_pcs"
    ),
    names(manifest)
  )
  if ("split" %in% names(manifest)) {
    manifest <- manifest[split %in% split_order]
  }
  manifest[, keep_cols, with = FALSE]
}

summarize_eqtl_analysis <- function(config, degs, g2sym, gwas_sets,
                                    repo_root, include_deg_views = FALSE,
                                    context_order = names(SEURAT_CONTEXT_TO_DATASET_ID),
                                    split_order = DEG_SPLITS) {
  require_data_table()
  gwas_sets <- normalize_gwas_sets(gwas_sets)
  gwas_disorders <- names(gwas_sets)
  manifest_qc <- manifest_qc_table(config$tqtl_in_dir, split_order = split_order)
  manifest <- load_eqtl_manifest(
    tqtl_in_dir = config$tqtl_in_dir,
    context_order = context_order,
    split_order = split_order
  )
  deg_global <- data.table::as.data.table(degs$global)

  summary_rows <- lapply(seq_len(nrow(manifest)), function(i) {
    row <- manifest[i]
    dataset_id <- row$dataset_id[[1]]
    context <- row$context[[1]]
    split <- row$split[[1]]

    map_cis_info <- read_map_cis(dataset_id, tqtl_out_dir = config$tqtl_out_dir, repo_root = repo_root)
    indep_info <- read_map_independent(dataset_id, tqtl_out_dir = config$tqtl_out_dir, repo_root = repo_root)

    list(
      map_cis = if (identical(map_cis_info$status, "ok")) {
        tag_eqtl_results(map_cis_info$dt, dataset_id, context, split, deg_global, g2sym, gwas_sets, repo_root = repo_root)
      } else {
        NULL
      },
      independent = if (identical(indep_info$status, "ok")) {
        tag_eqtl_results(indep_info$dt, dataset_id, context, split, deg_global, g2sym, gwas_sets, repo_root = repo_root)
      } else {
        NULL
      }
    )
  })

  map_cis_all <- data.table::rbindlist(lapply(summary_rows, `[[`, "map_cis"), use.names = TRUE, fill = TRUE)
  indep_all <- data.table::rbindlist(lapply(summary_rows, `[[`, "independent"), use.names = TRUE, fill = TRUE)

  map_cis_significant <- map_cis_all[qval < 0.05]

  map_cis_summary <- summarize_tagged_eqtls(
    map_cis_significant,
    context_order = context_order,
    split_order = split_order,
    gwas_disorders = gwas_disorders
  )

  parent_q <- map_cis_all[, .(phenotype_id, dataset_id, qval_parent = qval)]
  indep_f <- merge(indep_all, parent_q, by = c("phenotype_id", "dataset_id"), all = FALSE)
  indep_f <- indep_f[qval_parent < 0.05 & pval_perm < 0.05]
  indep_f <- order_annotated_eqtl_cols(indep_f)
  map_independent_summary <- summarize_independent_signals(
    indep_f,
    context_order = context_order,
    split_order = split_order,
    gwas_disorders = gwas_disorders
  )
  map_significant_unified <- unify_significant_eqtls(map_cis_significant, indep_f)
  map_significant_pairs <- collapse_significant_eqtl_pairs(map_significant_unified)
  map_significant_summary <- summarize_significant_pairs(
    map_significant_pairs,
    context_order = context_order,
    split_order = split_order,
    disorder_related = degs$disorder_related,
    gwas_disorders = gwas_disorders
  )

  deg_view_tables <- NULL
  if (isTRUE(include_deg_views)) {
    deg_view_tables <- list(
      map_cis_deg_view_summary = summarize_eqtl_deg_views(
        map_cis_significant,
        manifest = manifest,
        degs = degs,
        context_order = context_order,
        split_order = split_order,
        gwas_disorders = gwas_disorders
      ),
      map_independent_deg_view_summary = summarize_eqtl_deg_views(
        indep_f,
        manifest = manifest,
        degs = degs,
        context_order = context_order,
        split_order = split_order,
        signal_count = TRUE,
        gwas_disorders = gwas_disorders
      )
    )
  }

  list(
    analysis = config$analysis,
    manifest_qc = manifest_qc,
    manifest = manifest,
    map_cis_all = map_cis_all,
    indep_all = indep_all,
    indep_f = indep_f,
    map_cis_significant = map_cis_significant,
    map_independent_significant = indep_f,
    map_significant_unified = map_significant_unified,
    map_significant_pairs = map_significant_pairs,
    map_cis_summary = map_cis_summary,
    map_independent_summary = map_independent_summary,
    map_significant_summary = map_significant_summary,
    deg_view_tables = deg_view_tables
  )
}

assert_cols <- function(dt, cols, label = deparse(substitute(dt))) {
  missing_cols <- setdiff(cols, names(dt))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", label, ": ", paste(missing_cols, collapse = ", "))
  }
  invisible(TRUE)
}

first_non_missing <- function(x) {
  y <- x[!is.na(x) & nzchar(as.character(x))]
  if (length(y) == 0) return(NA_character_)
  as.character(y[[1]])
}

coloc_nominal_parquet_files <- function(dataset_id, tqtl_out_dir) {
  pattern <- paste0("^", gsub("\\.", "\\\\.", dataset_id), "\\.gene\\.cis_qtl_pairs\\.chr.*\\.parquet$")
  sort(list.files(tqtl_out_dir, pattern = pattern, full.names = TRUE))
}

read_coloc_nominal_dataset <- function(dataset_id, tqtl_out_dir, cis_window = 1000000L,
                                       chromosomes = NULL) {
  if (!requireNamespace("arrow", quietly = TRUE)) stop("The arrow package is required")
  if (!requireNamespace("tidyselect", quietly = TRUE)) stop("The tidyselect package is required")
  require_data_table()

  files <- coloc_nominal_parquet_files(dataset_id = dataset_id, tqtl_out_dir = tqtl_out_dir)
  if (!is.null(chromosomes)) {
    chr_pat <- paste0("\\.cis_qtl_pairs\\.(", paste(chromosomes, collapse = "|"), ")\\.parquet$")
    files <- files[grepl(chr_pat, files)]
  }
  if (length(files) == 0L) stop("No nominal parquet files found for dataset_id: ", dataset_id)

  cols <- c(
    "phenotype_id", "variant_id", "start_distance", "af", "ma_samples",
    "ma_count", "pval_nominal", "slope", "slope_se"
  )
  dt <- data.table::rbindlist(lapply(files, function(path) {
    as.data.table(arrow::read_parquet(path, col_select = tidyselect::all_of(cols)))
  }), use.names = TRUE, fill = TRUE)
  assert_cols(dt, cols, paste0(dataset_id, " nominal parquet"))
  dt[, `:=`(
    phenotype_id = as.character(phenotype_id),
    variant_id = as.character(variant_id),
    start_distance = as.integer(start_distance),
    af = as.numeric(af),
    pval_nominal = as.numeric(pval_nominal),
    slope = as.numeric(slope),
    slope_se = as.numeric(slope_se)
  )]
  dt <- dt[
    !is.na(phenotype_id) & nzchar(phenotype_id) &
      !is.na(variant_id) & nzchar(variant_id) &
      abs(start_distance) <= cis_window
  ]
  dt[]
}

coloc_gwas_dataset_cache_file <- function(dis, dataset_id, coloc_dir, si_min = 0.8) {
  dis <- gwas_check_disorder(dis)
  file.path(
    coloc_dir,
    dis,
    "cache",
    sprintf("gwas_%s_nominal-variants_SI%s.tsv.gz", dataset_id, gwas_si_tag(si_min))
  )
}

## convert nominal tensorQTL variant IDs into 1-bp BED regions.
## the PLINK2 pvar is the source of chromosome and position truth so GWAS
## queries use the same build and variant universe as tensorQTL.
coloc_regions_from_variants <- function(variant_ids, plink2_prefix) {
  require_data_table()
  info <- read_plink2_variant_table(plink2_prefix)
  out <- info[data.table::data.table(variant_id = unique(as.character(variant_ids))), on = "variant_id", nomatch = 0L]
  out <- out[!is.na(CHROM) & !is.na(POS), .(CHROM, start0 = POS - 1L, end1 = POS, variant_id)]
  data.table::setorder(out, CHROM, start0, end1, variant_id)
  out[]
}

## normalize raw bcftools output from the full GWAS BCF.
## ES/SE/LP/NE/NS/NC/SI become beta, beta_se, p, N, ns, ncas, and impinfo.
## this function applies the INFO score threshold only; coloc inputs must keep
## dense regional summary statistics without a GWAS p-value cutoff.
normalize_coloc_gwas_query_table <- function(dt, dis, si_min = 0.8) {
  dis <- gwas_check_disorder(dis)
  query_cols <- c("chr", "pos", "rsid", "a0", "a1", "beta", "beta_se", "lp", "N", "ns", "ncas", "impinfo")

  if (nrow(dt) == 0L) {
    dt <- data.table::as.data.table(stats::setNames(replicate(length(query_cols), logical(), simplify = FALSE), query_cols))
  } else {
    if (ncol(dt) != length(query_cols)) {
      stop("Unexpected coloc GWAS query column count: ", ncol(dt), "; expected ", length(query_cols))
    }
    data.table::setnames(dt, query_cols)
  }

  dt[, `:=`(
    chr = as.character(chr),
    pos = as.integer(pos),
    rsid = as.character(rsid),
    a0 = as.character(a0),
    a1 = as.character(a1),
    beta = suppressWarnings(as.numeric(beta)),
    beta_se = suppressWarnings(as.numeric(beta_se)),
    lp = suppressWarnings(as.numeric(lp)),
    N = suppressWarnings(as.numeric(N)),
    ns = suppressWarnings(as.numeric(ns)),
    ncas = suppressWarnings(as.numeric(ncas)),
    impinfo = suppressWarnings(as.numeric(impinfo))
  )]
  dt[, p := 10^(-lp)]
  dt[, variant_id := sprintf("%s:%s:%s:%s", chr, pos, a0, a1)]
  dt <- dt[is.na(impinfo) | impinfo >= si_min]
  data.table::setcolorder(dt, c(
    "rsid", "chr", "pos", "a0", "a1", "beta", "beta_se", "N",
    "ns", "ncas", "p", "impinfo", "lp", "variant_id"
  ))
  dt[]
}

## extract the dense GWAS slice needed for coloc for one disorder/domainCT.
## the function queries the disorder-specific full BCF at nominal eQTL variant
## positions, harmonizes alleles to PLINK2/tensorQTL IDs, flips beta for swapped
## alleles, writes a per-domainCT cache, and returns only matched variants.
extract_coloc_gwas_for_variants <- function(dis, variant_ids, plink2_prefix, out_file,
                                            repo_root = NULL, genotype_dir = NULL,
                                            si_min = 0.8, bcftools = "bcftools",
                                            use_cache = TRUE) {
  dis <- gwas_check_disorder(dis)
  if (isTRUE(use_cache) && file.exists(out_file)) return(data.table::fread(out_file))
  if (!nzchar(Sys.which(bcftools))) stop("Missing bcftools executable: ", bcftools)

  regions <- coloc_regions_from_variants(variant_ids = variant_ids, plink2_prefix = plink2_prefix)
  if (nrow(regions) == 0L) stop("No PLINK2 variants matched nominal variants for GWAS extraction.")

  bcf_file <- gwas_bcf_path(dis, repo_root = repo_root, genotype_dir = genotype_dir)
  region_file <- tempfile(pattern = "coloc-gwas-regions-", fileext = ".bed")
  query_file <- tempfile(pattern = "coloc-gwas-query-", fileext = ".tab")
  err_file <- tempfile(pattern = "coloc-gwas-query-err-", fileext = ".log")
  on.exit(unlink(c(region_file, query_file, err_file)), add = TRUE)

  ## query all nominal variant positions; do not impose a GWAS p-value cutoff.
  data.table::fwrite(unique(regions[, .(CHROM, start0, end1)]), region_file, sep = "\t", col.names = FALSE)
  status <- system2(
    bcftools,
    args = c(
      "query",
      "-R", shQuote(region_file),
      "-f", shQuote("%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%ES\t%SE\t%LP\t%NE\t%NS\t%NC\t%SI]\n"),
      shQuote(bcf_file)
    ),
    stdout = query_file,
    stderr = err_file
  )
  if (!identical(status, 0L)) {
    err <- if (file.exists(err_file)) readLines(err_file, warn = FALSE) else character()
    stop("bcftools dense GWAS query failed for ", dis, if (length(err)) paste0("\n", paste(err, collapse = "\n")) else "")
  }

  raw <- if (file.exists(query_file) && file.info(query_file)$size > 0) {
    data.table::fread(query_file, header = FALSE)
  } else {
    data.table::data.table()
  }
  gwas <- normalize_coloc_gwas_query_table(raw, dis = dis, si_min = si_min)
  matched <- match_coloc_gwas_geno(gwas, plink2_prefix = plink2_prefix)
  keep_cols <- c("variant_id", "rsid", "A1", "A2", "beta", "beta_se", "N", "p", "ncas", "impinfo", "match_mode")
  assert_cols(matched, keep_cols, paste0(dis, " dense matched GWAS"))
  matched <- matched[variant_id %in% unique(as.character(variant_ids))]

  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(matched, out_file, sep = "\t", quote = FALSE)
  matched[]
}

## align GWAS alleles to the PLINK2 pvar allele order used by tensorQTL.
## exact matches keep beta as-is; swapped REF/ALT matches flip beta so GWAS
## effects are on the same alternate allele encoded in the tensorQTL variant ID.
match_coloc_gwas_geno <- function(gwas, plink2_prefix) {
  require_data_table()
  req <- c("rsid", "chr", "pos", "a0", "a1", "beta", "beta_se", "N", "p", "ncas", "impinfo")
  assert_cols(gwas, req, "dense coloc GWAS")

  pvar_info <- read_plink2_pvar(plink2_prefix)
  gwas_dt <- data.table::copy(gwas)
  gwas_dt[, gwas_row_id := .I]

  exact <- merge(
    gwas_dt,
    pvar_info,
    by.x = c("chr", "pos", "a0", "a1"),
    by.y = c("chr", "pos", "pvar_ref", "pvar_alt"),
    all = FALSE,
    allow.cartesian = TRUE
  )
  if (nrow(exact) > 0L) {
    exact[, `:=`(
      variant_id = pvar_variant_id,
      A1 = a1,
      A2 = a0,
      match_mode = "exact",
      match_rank = 1L
    )]
  }

  swapped <- merge(
    gwas_dt,
    pvar_info,
    by.x = c("chr", "pos", "a0", "a1"),
    by.y = c("chr", "pos", "pvar_alt", "pvar_ref"),
    all = FALSE,
    allow.cartesian = TRUE
  )
  if (nrow(swapped) > 0L) {
    swapped[, `:=`(
      beta = -beta,
      variant_id = pvar_variant_id,
      A1 = a0,
      A2 = a1,
      match_mode = "swapped",
      match_rank = 2L
    )]
  }

  out_cols <- c("variant_id", "rsid", "A1", "A2", "beta", "beta_se", "N", "p", "ncas", "impinfo", "match_mode")
  out <- data.table::rbindlist(list(exact, swapped), use.names = TRUE, fill = TRUE)
  if (nrow(out) == 0L) {
    return(data.table::as.data.table(stats::setNames(replicate(length(out_cols), logical(), simplify = FALSE), out_cols)))
  }

  data.table::setorder(out, variant_id, match_rank, p)
  out[!duplicated(variant_id), ..out_cols]
}

## coloc case-control datasets need s, the case fraction.
## use the median ncas / N across matched GWAS rows to tolerate small row-level
## metadata variation in the source BCF.
coloc_case_fraction <- function(gwas_dt) {
  assert_cols(gwas_dt, c("N", "ncas"), "GWAS coloc table")
  vals <- gwas_dt[is.finite(N) & N > 0 & is.finite(ncas) & ncas > 0, ncas / N]
  vals <- vals[is.finite(vals) & vals > 0 & vals < 1]
  if (length(vals) == 0L) stop("Cannot infer GWAS case fraction from N and ncas")
  stats::median(vals, na.rm = TRUE)
}

coloc_assign_category_from_pp <- function(pp3, pp4) {
  pp34 <- pp3 + pp4
  ratio <- ifelse(!is.na(pp34) & pp34 > 0, pp4 / pp34, NA_real_)
  data.table::fcase(
    pp4 > 0.8, "strong_coloc",
    pp34 > 0.8 & ratio >= 0.9, "likely_coloc",
    pp34 > 0.8 & ratio < 0.5, "distinct_causal",
    (pp34 > 0.8) | (ratio > 0.8), "follow_up",
    default = NA_character_
  )
}

build_coloc_join_table <- function(eqtl_dt, gwas_dt, plink2_prefix) {
  require_data_table()
  assert_cols(eqtl_dt, c("phenotype_id", "variant_id", "start_distance", "af", "pval_nominal", "slope", "slope_se"), "nominal eQTL")
  assert_cols(gwas_dt, c("variant_id", "rsid", "beta", "beta_se", "N", "p", "ncas"), "matched GWAS")

  pvar <- read_plink2_variant_table(plink2_prefix)[, .(variant_id, CHROM, POS, REF, ALT)]
  e <- copy(eqtl_dt)
  e[, `:=`(
    snp = variant_id,
    beta_eqtl = as.numeric(slope),
    varbeta_eqtl = as.numeric(slope_se)^2,
    p_eqtl = as.numeric(pval_nominal),
    maf_eqtl = pmin(as.numeric(af), 1 - as.numeric(af))
  )]
  e <- merge(e, pvar, by = "variant_id", all.x = TRUE, sort = FALSE)

  g <- copy(gwas_dt)
  g[, `:=`(
    snp = variant_id,
    beta_gwas = as.numeric(beta),
    varbeta_gwas = as.numeric(beta_se)^2,
    p_gwas = as.numeric(p),
    n_gwas = as.numeric(N),
    ncas_gwas = as.numeric(ncas)
  )]

  out <- merge(
    e[, .(
      phenotype_id, snp, variant_id, pos = POS, start_distance,
      beta_eqtl, varbeta_eqtl, p_eqtl, maf_eqtl,
      ma_samples, ma_count, af
    )],
    g[, .(
      snp, rsid, beta_gwas, varbeta_gwas, p_gwas,
      n_gwas, ncas_gwas, impinfo, match_mode
    )],
    by = "snp",
    all = FALSE,
    sort = FALSE
  )
  out[!is.na(pos)]
}

run_coloc_abf_one_gene <- function(gene_id, coloc_dt, n_eqtl, s_gwas,
                                   p1 = 1e-4, p2 = 1e-4, p12 = 1e-5,
                                   min_snps = 10L, min_abs_eqtl_z = 2,
                                   warn_minp = 1e-6) {
  if (!requireNamespace("coloc", quietly = TRUE)) stop("The coloc package is required")
  gene_id_value <- gene_id
  d <- coloc_dt[phenotype_id == gene_id_value]
  d <- d[
    is.finite(beta_eqtl) & is.finite(varbeta_eqtl) & varbeta_eqtl > 0 &
      is.finite(beta_gwas) & is.finite(varbeta_gwas) & varbeta_gwas > 0 &
      is.finite(p_eqtl) & p_eqtl > 0 & p_eqtl <= 1 &
      is.finite(p_gwas) & p_gwas > 0 & p_gwas <= 1 &
      is.finite(maf_eqtl) & maf_eqtl > 0 & maf_eqtl < 1
  ]
  if (anyDuplicated(d$snp)) {
    data.table::setorder(d, snp, p_eqtl, p_gwas)
    d <- d[!duplicated(snp)]
  }
  if (nrow(d) < min_snps) return(NULL)

  z_eqtl <- abs(d$beta_eqtl / sqrt(d$varbeta_eqtl))
  if (!any(is.finite(z_eqtl) & z_eqtl >= min_abs_eqtl_z)) return(NULL)

  eqtl <- list(
    snp = d$snp,
    position = d$pos,
    beta = d$beta_eqtl,
    varbeta = d$varbeta_eqtl,
    pvalues = d$p_eqtl,
    MAF = d$maf_eqtl,
    N = n_eqtl,
    type = "quant"
  )
  gwas <- list(
    snp = d$snp,
    position = d$pos,
    beta = d$beta_gwas,
    varbeta = d$varbeta_gwas,
    pvalues = d$p_gwas,
    N = stats::median(d$n_gwas, na.rm = TRUE),
    type = "cc",
    s = s_gwas
  )

  coloc::check_dataset(eqtl, warn.minp = warn_minp)
  coloc::check_dataset(gwas, warn.minp = warn_minp)
  res <- NULL
  invisible(utils::capture.output(
    res <- coloc::coloc.abf(dataset1 = eqtl, dataset2 = gwas, p1 = p1, p2 = p2, p12 = p12)
  ))
  res
}

run_coloc_abf_one_gene_safe <- function(gene_id, coloc_dt, n_eqtl, s_gwas, ...) {
  tryCatch(
    run_coloc_abf_one_gene(gene_id = gene_id, coloc_dt = coloc_dt, n_eqtl = n_eqtl, s_gwas = s_gwas, ...),
    error = function(e) structure(list(.error = conditionMessage(e)), class = "coloc_err")
  )
}

run_coloc_sensitivity_one <- function(obj, rule, npoints = 100L) {
  if (!requireNamespace("coloc", quietly = TRUE)) stop("The coloc package is required")
  tryCatch(
    {
      sens <- NULL
      invisible(utils::capture.output(
        withCallingHandlers(
          sens <- coloc::sensitivity(obj, rule = rule, doplot = FALSE, npoints = npoints),
          message = function(m) invokeRestart("muffleMessage")
        )
      ))
      as.data.table(sens)
    },
    error = function(e) data.table::data.table(error = conditionMessage(e))
  )
}

summarize_coloc_sensitivity <- function(obj, rule = "H4 > 0.8", npoints = 100L,
                                        p12_plausible_min = 5e-6,
                                        p12_plausible_max = 5e-5) {
  sens <- run_coloc_sensitivity_one(obj, rule = rule, npoints = npoints)
  if ("error" %in% names(sens)) {
    return(data.table::data.table(
      sensitivity_error = sens$error[[1]],
      npoints = NA_integer_,
      n_pass = NA_integer_,
      frac_pass = NA_real_,
      n_plausible = NA_integer_,
      n_plausible_pass = NA_integer_,
      frac_plausible_pass = NA_real_,
      min_pass_p12 = NA_real_,
      max_pass_p12 = NA_real_,
      rule = rule
    ))
  }
  pass_p12 <- sens[pass == TRUE, p12]
  plausible <- sens[p12 >= p12_plausible_min & p12 <= p12_plausible_max]
  data.table::data.table(
    sensitivity_error = NA_character_,
    npoints = nrow(sens),
    n_pass = sens[, sum(pass == TRUE, na.rm = TRUE)],
    frac_pass = sens[, mean(pass == TRUE, na.rm = TRUE)],
    n_plausible = nrow(plausible),
    n_plausible_pass = plausible[, sum(pass == TRUE, na.rm = TRUE)],
    frac_plausible_pass = if (nrow(plausible) > 0L) plausible[, mean(pass == TRUE, na.rm = TRUE)] else NA_real_,
    min_pass_p12 = if (length(pass_p12)) min(pass_p12, na.rm = TRUE) else NA_real_,
    max_pass_p12 = if (length(pass_p12)) max(pass_p12, na.rm = TRUE) else NA_real_,
    rule = rule
  )
}

flatten_coloc_results <- function(result_list, disorder, dataset_id, context,
                                  split = "all", g2sym = NULL, cs_level = 0.95,
                                  candidate_only = FALSE) {
  require_data_table()
  empty_out <- function() data.table::data.table(
    disorder = character(),
    dataset_id = character(),
    context = character(),
    split = character(),
    gene_id = character(),
    gene_name = character(),
    cat = character(),
    nsnps = numeric(),
    PP0 = numeric(),
    PP1 = numeric(),
    PP2 = numeric(),
    PP3 = numeric(),
    PP4 = numeric(),
    PP34 = numeric(),
    PP4_over_PP34 = numeric(),
    lead_snp = character(),
    lead_snp_PPH4 = numeric(),
    lead_snp_PPsho = numeric(),
    cs95_n_snp = integer(),
    p1 = numeric(),
    p2 = numeric(),
    p12 = numeric()
  )
  if (length(result_list) == 0L) return(empty_out())

  disorder_value <- disorder
  dataset_id_value <- dataset_id
  context_value <- context
  split_value <- split
  rows <- lapply(names(result_list), function(gene_id) {
    gene_id_value <- gene_id
    obj <- result_list[[gene_id]]
    if (is.null(obj) || is.null(obj$summary) || is.null(obj$results)) return(NULL)
    sm <- as.list(obj$summary)
    res <- data.table::as.data.table(obj$results)
    if (!"SNP.PP.H4" %in% names(res)) return(NULL)
    data.table::setorderv(res, "SNP.PP.H4", -1L)
    res[, cum_h4 := cumsum(SNP.PP.H4)]
    cs_idx <- which(res$cum_h4 >= cs_level)[1L]
    cs_n <- if (is.na(cs_idx)) nrow(res) else cs_idx
    lead <- res[1]
    pp3 <- as.numeric(sm[["PP.H3.abf"]])
    pp4 <- as.numeric(sm[["PP.H4.abf"]])
    pp34 <- pp3 + pp4
    ratio <- if (!is.na(pp34) && pp34 > 0) pp4 / pp34 else NA_real_

    data.table::data.table(
      disorder = disorder_value,
      dataset_id = dataset_id_value,
      context = context_value,
      split = split_value,
      gene_id = gene_id_value,
      gene_name = if (!is.null(g2sym)) unname(g2sym[[gene_id_value]]) else NA_character_,
      cat = coloc_assign_category_from_pp(pp3, pp4),
      nsnps = as.numeric(sm[["nsnps"]]),
      PP0 = as.numeric(sm[["PP.H0.abf"]]),
      PP1 = as.numeric(sm[["PP.H1.abf"]]),
      PP2 = as.numeric(sm[["PP.H2.abf"]]),
      PP3 = pp3,
      PP4 = pp4,
      PP34 = pp34,
      PP4_over_PP34 = ratio,
      lead_snp = as.character(lead[["snp"]]),
      lead_snp_PPH4 = as.numeric(lead[["SNP.PP.H4"]]),
      lead_snp_PPsho = pp4 * as.numeric(lead[["SNP.PP.H4"]]),
      cs95_n_snp = as.integer(cs_n),
      p1 = as.numeric(obj$priors[["p1"]]),
      p2 = as.numeric(obj$priors[["p2"]]),
      p12 = as.numeric(obj$priors[["p12"]])
    )
  })
  out <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (is.null(out) || ncol(out) == 0L) return(empty_out())
  if (isTRUE(candidate_only) && nrow(out) > 0L) out <- out[!is.na(cat)]
  out[]
}

flatten_coloc_snps <- function(result_list, disorder, dataset_id, context, split = "all") {
  require_data_table()
  disorder_value <- disorder
  dataset_id_value <- dataset_id
  context_value <- context
  split_value <- split
  rows <- lapply(names(result_list), function(gene_id) {
    gene_id_value <- gene_id
    obj <- result_list[[gene_id]]
    if (is.null(obj) || is.null(obj$results)) return(NULL)
    dt <- data.table::as.data.table(obj$results)
    dt[, `:=`(
      disorder = disorder_value,
      dataset_id = dataset_id_value,
      context = context_value,
      split = split_value,
      gene_id = gene_id_value
    )]
    data.table::setcolorder(dt, c(
      "disorder", "dataset_id", "context", "split", "gene_id",
      setdiff(names(dt), c("disorder", "dataset_id", "context", "split", "gene_id"))
    ))
    dt
  })
  data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}
