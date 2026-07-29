library(data.table)
library(testthat)

utils_path <- c("utils.R", "../utils.R", "../../utils.R")
utils_path <- utils_path[file.exists(utils_path)][1]
source(utils_path, chdir = FALSE)

with_gwas_attrs <- function(dt, cache_dir = tempdir()) {
  attr(dt, "gwas_dis") <- "SCZD"
  attr(dt, "gwas_pval") <- 1e-6
  attr(dt, "gwas_si_min") <- 0.8
  attr(dt, "gwas_cache_dir") <- cache_dir
  dt
}

test_that("GWAS cache names are stable", {
  genotype_dir <- tempdir()

  expect_equal(gwas_pval_tag(1e-6), "1e-6")
  expect_equal(gwas_pval_tag(5e-8), "5e-8")

  cache_file <- gwas_cache_file("SCZD", 1e-6, genotype_dir = genotype_dir)
  expect_equal(
    basename(cache_file),
    "GWAS-SCZD_flt_p1e-6_SI0.8.hg38.tab.gz"
  )
  expect_equal(
    basename(gwas_cache_file("MDD", 1e-5, genotype_dir = genotype_dir)),
    "GWAS-MDD_full23andMe_flt_p1e-5_SInone.hg38.tab.gz"
  )
  expect_true(is.na(gwas_resolve_si_min("MDD")))
  expect_equal(gwas_resolve_si_min("SCZD"), 0.8)
})

test_that("disorder names resolve to BCF paths", {
  genotype_dir <- file.path(tempdir(), "geno-bcf")
  bcf_file <- file.path(genotype_dir, GWAS_BCF_FILES[["SCZD"]])
  dir.create(dirname(bcf_file), recursive = TRUE, showWarnings = FALSE)
  file.create(bcf_file)

  expect_equal(gwas_bcf_path("SCZD", genotype_dir = genotype_dir), normalizePath(bcf_file))
  expect_error(gwas_bcf_path("BAD", genotype_dir = genotype_dir), "Unsupported disorder")
})

test_that("GWAS query output is normalized", {
  raw <- fread(
    text = "chr1\t100\trs1\tA\tG\t0.2\t0.1\t6\t1000\t2000\t500\t0.9",
    header = FALSE
  )

  gwas <- normalize_gwas_query_table(raw, dis = "SCZD", pval = 1e-6, si_min = 0.8)

  expect_equal(gwas$rsid, "rs1")
  expect_equal(gwas$variant_id, "chr1:100:A:G")
  expect_equal(gwas$p, 1e-6, tolerance = 1e-12)
  expect_equal(gwas$impinfo, 0.9)
})

test_that("loadGWAS reads an existing cache without bcftools", {
  genotype_dir <- file.path(tempdir(), "geno-cache")
  dir.create(genotype_dir, recursive = TRUE, showWarnings = FALSE)
  cache_file <- gwas_cache_file("SCZD", 1e-6, genotype_dir = genotype_dir)
  fwrite(
    data.table(
      rsid = "rs-cache", chr = "chr1", pos = 1L, a0 = "A", a1 = "G",
      beta = 0.1, beta_se = 0.01, N = 100, p = 1e-7,
      impinfo = 0.9, ncas = 50, ns = 100, lp = 7,
      variant_id = "chr1:1:A:G"
    ),
    cache_file,
    sep = "\t"
  )

  gwas <- loadGWAS("SCZD", 1e-6, genotype_dir = genotype_dir, bcftools = "missing-bcftools")

  expect_equal(gwas$rsid, "rs-cache")
  expect_equal(attr(gwas, "gwas_dis"), "SCZD")
})

test_that("loadGWAS can use bcftools query output", {
  genotype_dir <- file.path(tempdir(), "geno-query")
  bcf_file <- file.path(genotype_dir, GWAS_BCF_FILES[["SCZD"]])
  dir.create(dirname(bcf_file), recursive = TRUE, showWarnings = FALSE)
  file.create(bcf_file)

  fake_bcftools <- file.path(tempdir(), "fake-bcftools")
  writeLines(
    c(
      "#!/bin/sh",
      "printf '%s\\n' 'chr1\t100\trs1\tA\tG\t0.2\t0.1\t6\t1000\t2000\t500\t0.9'"
    ),
    fake_bcftools
  )
  Sys.chmod(fake_bcftools, "0755")

  gwas <- loadGWAS(
    "SCZD", 1e-6,
    genotype_dir = genotype_dir,
    bcftools = fake_bcftools,
    use_cache = FALSE
  )

  expect_equal(gwas$variant_id, "chr1:100:A:G")
  expect_true(file.exists(gwas_cache_file("SCZD", 1e-6, genotype_dir = genotype_dir)))
})

test_that("matchGwasGeno handles exact, swapped, and unmatched variants", {
  cache_dir <- file.path(tempdir(), "geno-match")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  plink_prefix <- file.path(cache_dir, "merged_test")
  writeLines(
    c(
      "##filedate=20260506",
      "#CHROM\tPOS\tID\tREF\tALT\tINFO",
      "chr1\t100\tchr1:100:A:G\tA\tG\t.",
      "chr1\t200\tchr1:200:T:C\tT\tC\t.",
      "chr1\t300\tchr1:300:A:C\tA\tC\t."
    ),
    paste0(plink_prefix, ".pvar")
  )

  gwas <- with_gwas_attrs(
    data.table(
      rsid = c("rs-exact", "rs-swap", "rs-drop"),
      chr = c("chr1", "chr1", "chr1"),
      pos = c(100L, 200L, 400L),
      a0 = c("A", "C", "A"),
      a1 = c("G", "T", "G"),
      beta = c(0.2, 0.4, 0.6),
      beta_se = c(0.01, 0.02, 0.03),
      N = c(100, 100, 100),
      p = c(1e-7, 1e-8, 1e-9),
      ncas = c(50, 50, 50),
      impinfo = c(0.9, 0.9, 0.9)
    ),
    cache_dir = cache_dir
  )

  matched <- matchGwasGeno(gwas, plink_prefix, use_cache = FALSE)

  expect_equal(nrow(matched), 2L)
  expect_equal(matched[rsid == "rs-exact", variant_id], "chr1:100:A:G")
  expect_equal(matched[rsid == "rs-exact", beta], 0.2)
  expect_equal(matched[rsid == "rs-exact", match_mode], "exact")
  expect_equal(matched[rsid == "rs-swap", variant_id], "chr1:200:T:C")
  expect_equal(matched[rsid == "rs-swap", beta], -0.4)
  expect_equal(matched[rsid == "rs-swap", A1], "C")
  expect_false("rs-drop" %in% matched$rsid)
})

test_that("matchGwasGeno reads an existing cache before pvar", {
  cache_dir <- file.path(tempdir(), "geno-match-cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  gwas <- with_gwas_attrs(data.table(rsid = "rs1"), cache_dir = cache_dir)
  plink_prefix <- file.path(cache_dir, "missing_pvar")
  cache_file <- gwas_tqtl_cache_file(gwas, plink_prefix)
  fwrite(
    data.table(
      variant_id = "chr1:1:A:G", rsid = "rs1", A1 = "G", A2 = "A",
      beta = 0.1, beta_se = 0.01, N = 100, p = 1e-7,
      ncas = 50, impinfo = 0.9, match_mode = "exact"
    ),
    cache_file,
    sep = "\t"
  )

  matched <- matchGwasGeno(gwas, plink_prefix)

  expect_equal(matched$variant_id, "chr1:1:A:G")
})
