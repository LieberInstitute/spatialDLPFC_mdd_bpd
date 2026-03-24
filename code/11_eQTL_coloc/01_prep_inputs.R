#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(edgeR)
  library(matrixStats)
  library(sva)
  library(qs2)
})

get_script_path <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", ca, value = TRUE)
  if (length(m) == 0) stop("Cannot determine script path from commandArgs")
  sub("^--file=", "", m[[1]])
}

script_dir <- dirname(normalizePath(get_script_path()))
repo_root <- normalizePath(file.path(script_dir, "../.."))

opt <- list(
  spe_file = file.path(repo_root, "processed-data", "06_pseudobulk", "Seurat",
                       "spe_n119_pseudo-no-lowUMI_sample-seurat-pc30_norm-filt.Rdata"),
  geno_prefix = file.path(repo_root, "processed-data", "00_genotypes", "plink2", "merged_maf05"),
  snp_pcs_file = file.path(repo_root, "processed-data", "00_genotypes", "plink2", "merged_maf05_pca.eigenvec"),
  out_dir = file.path(repo_root, "processed-data", "11_eQTL_coloc", "tqtl_in"),
  splits = "all,male,female",
  min_samples = 20L,
  write_spe = FALSE,
  check_only = FALSE
)

args <- commandArgs(trailingOnly = TRUE)
i <- 1L
while (i <= length(args)) {
  a <- args[[i]]
  if (!startsWith(a, "--")) stop("Unexpected argument: ", a)
  key <- sub("^--", "", a)
  if (key %in% c("write-spe", "check-only")) {
    if (key == "write-spe") opt$write_spe <- TRUE
    if (key == "check-only") opt$check_only <- TRUE
    i <- i + 1L
    next
  }
  if (i == length(args)) stop("Missing value for argument: ", a)
  val <- args[[i + 1L]]
  if (key == "spe-file") opt$spe_file <- val
  else if (key == "geno-prefix") opt$geno_prefix <- val
  else if (key == "snp-pcs-file") opt$snp_pcs_file <- val
  else if (key == "out-dir") opt$out_dir <- val
  else if (key == "splits") opt$splits <- val
  else if (key == "min-samples") opt$min_samples <- as.integer(val)
  else stop("Unknown argument: --", key)
  i <- i + 2L
}

if (!file.exists(opt$spe_file)) stop("Missing spe file: ", opt$spe_file)
if (!file.exists(paste0(opt$geno_prefix, ".pgen"))) stop("Missing genotype pgen: ", opt$geno_prefix, ".pgen")
if (!file.exists(paste0(opt$geno_prefix, ".psam"))) stop("Missing genotype psam: ", opt$geno_prefix, ".psam")
if (!file.exists(paste0(opt$geno_prefix, ".pvar"))) stop("Missing genotype pvar: ", opt$geno_prefix, ".pvar")
if (!file.exists(opt$snp_pcs_file)) stop("Missing SNP PCs file: ", opt$snp_pcs_file)
if (!dir.exists(opt$out_dir) && !opt$check_only) dir.create(opt$out_dir, recursive = TRUE)

cat("Repo root:", repo_root, "\n")
cat("Input SPE:", opt$spe_file, "\n")
cat("SNP PCs:", opt$snp_pcs_file, "\n")
cat("Output dir:", opt$out_dir, "\n")

spe_env <- new.env(parent = emptyenv())
load(opt$spe_file, envir = spe_env)
obj_names <- ls(spe_env)
if (length(obj_names) == 0) stop("No objects found in ", opt$spe_file)

is_spe <- vapply(obj_names, function(x) inherits(get(x, envir = spe_env), "SpatialExperiment"), logical(1))
if (!any(is_spe)) stop("No SpatialExperiment object found in ", opt$spe_file)
spe <- get(obj_names[which(is_spe)[1]], envir = spe_env)

cd <- as.data.frame(colData(spe))

if (!"seurat_label" %in% names(cd)) {
  seurat_cols <- grep("^seurat_", names(cd), value = TRUE)
  if (length(seurat_cols) == 0) stop("Missing seurat_label and no seurat_* column found")
  cd$seurat_label <- cd[[seurat_cols[[1]]]
  ]
}
cd$seurat_label <- as.factor(as.character(cd$seurat_label))

if (!"brnum" %in% names(cd)) {
  if ("BrNum" %in% names(cd)) {
    cd$brnum <- cd$BrNum
  } else {
    cd$brnum <- colnames(spe)
  }
}
cd$brnum <- as.character(cd$brnum)

if ("condition" %in% names(cd)) {
  dx0 <- toupper(as.character(cd$condition))
} else if ("DX" %in% names(cd)) {
  dx0 <- toupper(as.character(cd$DX))
} else if ("dx" %in% names(cd)) {
  dx0 <- toupper(as.character(cd$dx))
} else {
  stop("Missing diagnosis column: expected one of condition/DX/dx")
}
cd$DX <- factor(dx0, levels = c("NTC", "MDD", "BPD"))
if (any(is.na(cd$DX))) {
  bad <- unique(dx0[is.na(cd$DX)])
  stop("Unexpected diagnosis labels: ", paste(bad, collapse = ", "))
}

if (!"sex" %in% names(cd)) stop("Missing sex column")
cd$sex <- factor(toupper(as.character(cd$sex)), levels = c("M", "F"))
if (any(is.na(cd$sex))) stop("Sex column contains non M/F values")

if (!"age" %in% names(cd)) stop("Missing age column")
cd$age <- suppressWarnings(as.numeric(cd$age))
if (any(is.na(cd$age))) stop("age column contains non-numeric values")

if ("PC3" %in% names(cd)) {
  cd$PC3 <- suppressWarnings(as.numeric(cd$PC3))
} else if ("pc3" %in% names(cd)) {
  cd$PC3 <- suppressWarnings(as.numeric(cd$pc3))
} else if ("PCA_1663" %in% reducedDimNames(spe) && "PC3" %in% colnames(reducedDim(spe, "PCA_1663"))) {
  cd$PC3 <- reducedDim(spe, "PCA_1663")[, "PC3"]
} else if (length(reducedDimNames(spe)) > 0 && ncol(reducedDim(spe, reducedDimNames(spe)[1])) >= 3) {
  cd$PC3 <- reducedDim(spe, reducedDimNames(spe)[1])[, 3]
} else {
  stop("Could not derive PC3 from colData or reducedDims")
}
if (any(is.na(cd$PC3))) stop("PC3 contains NA values")

colData(spe)$seurat_label <- cd$seurat_label
colData(spe)$brnum <- cd$brnum
colData(spe)$DX <- cd$DX
colData(spe)$sex <- cd$sex
colData(spe)$age <- cd$age
colData(spe)$PC3 <- cd$PC3

snpPCs <- fread(opt$snp_pcs_file, data.table = FALSE)
if (!("SAMPLE_ID" %in% colnames(snpPCs))) {
  nidcols <- ncol(snpPCs) - length(grep("PC\\d+$", colnames(snpPCs)))
  if (nidcols == 2) snpPCs <- snpPCs[, -1, drop = FALSE]
  colnames(snpPCs)[1] <- "SAMPLE_ID"
}
rownames(snpPCs) <- snpPCs$SAMPLE_ID
colnames(snpPCs) <- gsub("^PC", "snpPC", colnames(snpPCs))
for (pc in paste0("snpPC", 1:5)) {
  if (!(pc %in% colnames(snpPCs))) stop("Missing ", pc, " in SNP PCs file")
}

prep_spe_qtl <- function(
  spe,
  bed_detect_prop = 0.20,
  bed_min_cpm = 0.10,
  bed_min_count = 6,
  min_detect_n_floor = 3,
  pca_pi0_max = 0.90,
  pca_use_hvg = TRUE,
  pca_use_detect_filter = TRUE,
  pca_hvg_n = 5000,
  prior_count = 1
) {
  stopifnot("counts" %in% assayNames(spe))

  cnt <- assay(spe, "counts")
  dge0 <- DGEList(cnt)
  dge0 <- calcNormFactors(dge0, method = "TMM")
  cpm0 <- cpm(dge0, log = FALSE)

  n <- ncol(spe)
  min_detect_n <- max(ceiling(bed_detect_prop * n), min_detect_n_floor)

  det <- (cpm0 >= bed_min_cpm) & (cnt >= bed_min_count)
  det_n <- rowSums(det)

  keep_gene_bed <- det_n >= min_detect_n
  spe_bed <- spe[keep_gene_bed, , drop = FALSE]

  dge_bed <- DGEList(assay(spe_bed, "counts"))
  dge_bed <- calcNormFactors(dge_bed, method = "TMM")
  log2cpm <- cpm(dge_bed, log = TRUE, prior.count = prior_count)

  rr <- rowRanks(log2cpm, ties.method = "average")
  pp <- (rr - 0.5) / ncol(log2cpm)
  log2cpm_rint <- qnorm(pp)

  assays(spe_bed)$tmm <- log2cpm
  assays(spe_bed)$rint <- log2cpm_rint

  cnt_qc <- assay(spe, "counts")
  pi0 <- rowMeans(cnt_qc == 0)
  keep_gene_pca <- pi0 < pca_pi0_max
  if (pca_use_detect_filter) keep_gene_pca <- keep_gene_pca & (det_n >= min_detect_n)

  spe_pca <- spe[keep_gene_pca, , drop = FALSE]

  dge_pca <- DGEList(assay(spe_pca, "counts"))
  dge_pca <- calcNormFactors(dge_pca, method = "TMM")
  log2cpm_pca <- cpm(dge_pca, log = TRUE, prior.count = prior_count)

  if (pca_use_hvg) {
    cpm_lin <- cpm(dge_pca, log = FALSE)
    mu <- rowMeans(cpm_lin)
    va <- rowVars(cpm_lin)
    fano <- va / pmax(mu, 1e-8)
    o <- order(fano, decreasing = TRUE)
    sel <- o[seq_len(min(pca_hvg_n, length(o)))]
    log2cpm_pca <- log2cpm_pca[sel, , drop = FALSE]
    spe_pca <- spe_pca[sel, , drop = FALSE]
  }

  mat_z <- t(scale(t(log2cpm_pca), center = TRUE, scale = TRUE))
  mat_z[!is.finite(mat_z)] <- 0

  assays(spe_pca)$tmm <- log2cpm_pca
  assays(spe_pca)$zscore <- mat_z

  list(
    spe_bed = spe_bed,
    spe_pca = spe_pca,
    n_genes_bed = nrow(spe_bed),
    n_genes_pca = nrow(spe_pca)
  )
}

spe2bed <- function(spe_bed) {
  rr_df <- as.data.frame(rowRanges(spe_bed))
  tmm <- assays(spe_bed)$tmm
  stopifnot(!is.null(tmm))

  rinvnorm <- function(x) qnorm((rank(x, ties.method = "average") - 0.5) / length(x))
  counts <- t(apply(tmm, 1, rinvnorm))

  rr_df$tss_start <- ifelse(rr_df$strand == "-", rr_df$end, rr_df$start)
  rr_df$start <- rr_df$tss_start
  rr_df <- rr_df %>%
    tibble::rownames_to_column("ID") %>%
    dplyr::arrange(seqnames, start) %>%
    dplyr::mutate(end = start + 1) %>%
    dplyr::select(`#Chr` = seqnames, start, end, ID)

  colnames(counts) <- colnames(spe_bed)
  as.data.frame(counts) %>%
    tibble::rownames_to_column("ID") %>%
    dplyr::left_join(rr_df, ., by = "ID")
}

covar_format <- function(data) {
  data <- as.data.frame(data)
  data <- t(data)
  as.data.frame(data) %>% rownames_to_column("id")
}

sanitize <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "-", x)
  x <- gsub("(^-+|-+$)", "", x)
  tolower(x)
}

dataset_id_from <- function(cluster, split) {
  base <- sanitize(cluster)
  if (base == "micro-vasc") base <- "uvasc"
  suffix <- switch(
    split,
    all = "",
    male = "_m",
    female = "_f",
    stop("Unsupported split: ", split)
  )
  paste0(base, suffix)
}

splits <- trimws(unlist(strsplit(opt$splits, ",", fixed = TRUE)))
splits <- splits[nchar(splits) > 0]
allowed_splits <- c("all", "male", "female")
if (!all(splits %in% allowed_splits)) stop("Unsupported split(s): ", paste(setdiff(splits, allowed_splits), collapse = ", "))

clusters <- sort(unique(as.character(colData(spe)$seurat_label)))
manifest <- list()

for (cluster in clusters) {
  cluster_mask <- as.character(colData(spe)$seurat_label) == cluster
  for (split in splits) {
    split_mask <- rep(TRUE, ncol(spe))
    if (split == "male") split_mask <- as.character(colData(spe)$sex) == "M"
    if (split == "female") split_mask <- as.character(colData(spe)$sex) == "F"

    sel <- cluster_mask & split_mask
    ds_id <- dataset_id_from(cluster, split)

    if (sum(sel) < opt$min_samples) {
      manifest[[length(manifest) + 1L]] <- data.frame(
        dataset_id = ds_id,
        seurat_label = cluster,
        split = split,
        n_samples = sum(sel),
        n_genes_bed = NA_integer_,
        n_genes_pca = NA_integer_,
        n_expr_pcs = NA_integer_,
        status = "skipped",
        reason = paste0("n_samples<", opt$min_samples),
        stringsAsFactors = FALSE
      )
      next
    }

    spe_sub <- spe[, sel, drop = FALSE]

    br <- as.character(colData(spe_sub)$brnum)
    if (any(is.na(br) | br == "")) {
      manifest[[length(manifest) + 1L]] <- data.frame(
        dataset_id = ds_id,
        seurat_label = cluster,
        split = split,
        n_samples = ncol(spe_sub),
        n_genes_bed = NA_integer_,
        n_genes_pca = NA_integer_,
        n_expr_pcs = NA_integer_,
        status = "skipped",
        reason = "missing brnum",
        stringsAsFactors = FALSE
      )
      next
    }

    colnames(spe_sub) <- br
    if (anyDuplicated(colnames(spe_sub))) {
      manifest[[length(manifest) + 1L]] <- data.frame(
        dataset_id = ds_id,
        seurat_label = cluster,
        split = split,
        n_samples = ncol(spe_sub),
        n_genes_bed = NA_integer_,
        n_genes_pca = NA_integer_,
        n_expr_pcs = NA_integer_,
        status = "skipped",
        reason = "duplicate sample IDs after brnum assignment",
        stringsAsFactors = FALSE
      )
      next
    }

    if (!all(colnames(spe_sub) %in% rownames(snpPCs))) {
      missing_ids <- setdiff(colnames(spe_sub), rownames(snpPCs))
      stop("Missing SNP PCs for sample IDs: ", paste(missing_ids, collapse = ", "))
    }

    for (pc in paste0("snpPC", 1:5)) {
      colData(spe_sub)[[pc]] <- as.numeric(snpPCs[colnames(spe_sub), pc])
    }

    pd <- as.data.frame(colData(spe_sub))
    pd$DX <- factor(as.character(pd$DX), levels = c("NTC", "MDD", "BPD"))
    pd$age <- as.numeric(pd$age)
    pd$PC3 <- as.numeric(pd$PC3)

    if (split == "all") {
      pd$sex <- factor(as.character(pd$sex), levels = c("M", "F"))
      model_formula <- as.formula("~ DX + sex + age + PC3 + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5")
    } else {
      model_formula <- as.formula("~ DX + age + PC3 + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5")
    }

    model <- model.matrix(model_formula, data = pd)

    rr_chr <- as.character(seqnames(rowRanges(spe_sub)))
    keep_chr <- !grepl("^(chr)?(Y|M|MT)$", rr_chr, ignore.case = TRUE)
    spe_sub <- spe_sub[keep_chr, , drop = FALSE]

    res <- prep_spe_qtl(spe_sub)

    lmx <- assays(res$spe_pca)$zscore
    n_expr_pcs <- 0L
    ffPCs <- matrix(numeric(), nrow = ncol(spe_sub), ncol = 0,
      dimnames = list(colnames(spe_sub), character()))

    if (nrow(lmx) > 1 && ncol(lmx) > 1) {
      pca <- prcomp(t(lmx))
      k <- tryCatch(as.integer(num.sv(lmx, model)), error = function(e) 0L)
      k <- max(0L, min(k, ncol(pca$x)))
      if (k > 0L) {
        ffPCs <- pca$x[, seq_len(k), drop = FALSE]
        colnames(ffPCs) <- paste0("exprPC", seq_len(ncol(ffPCs)))
        rownames(ffPCs) <- colnames(spe_sub)
        n_expr_pcs <- ncol(ffPCs)
      }
    }

    if (!opt$check_only) {
      fn_base <- file.path(opt$out_dir, paste0(ds_id, ".gene"))
      fn_exprpcs <- paste0(fn_base, ".exprPCs.qs2")
      fn_covars <- paste0(fn_base, ".covars.txt")
      fn_bed <- paste0(fn_base, ".expr.bed.gz")
      fn_spe <- paste0(fn_base, ".spe.qs2")

      qs_save(ffPCs, file = fn_exprpcs)

      cpd <- covar_format(model)
      if (ncol(ffPCs) > 0) {
        cpc <- covar_format(ffPCs)
        stopifnot(identical(colnames(cpd), colnames(cpc)))
        covars <- rbind(cpd, cpc)
      } else {
        covars <- cpd
      }
      fwrite(covars, file = fn_covars, sep = "\t", quote = FALSE, row.names = FALSE)

      bed <- spe2bed(res$spe_bed)
      fwrite(bed, file = fn_bed, sep = "\t", quote = FALSE, row.names = FALSE)

      if (opt$write_spe) qs_save(res$spe_bed, file = fn_spe)
    }

    manifest[[length(manifest) + 1L]] <- data.frame(
      dataset_id = ds_id,
      seurat_label = cluster,
      split = split,
      n_samples = ncol(spe_sub),
      n_genes_bed = res$n_genes_bed,
      n_genes_pca = res$n_genes_pca,
      n_expr_pcs = n_expr_pcs,
      status = ifelse(opt$check_only, "checked", "prepared"),
      reason = "",
      stringsAsFactors = FALSE
    )

    cat(sprintf("[%s] %s (%s): n=%d, genes_bed=%d, exprPCs=%d\n",
      ifelse(opt$check_only, "CHECK", "DONE"), ds_id, cluster, ncol(spe_sub), res$n_genes_bed, n_expr_pcs))
  }
}

manifest_df <- dplyr::bind_rows(manifest)
if (nrow(manifest_df) == 0) stop("No datasets were processed")

if (!opt$check_only) {
  manifest_file <- file.path(opt$out_dir, "prep_manifest.csv")
  fwrite(manifest_df, file = manifest_file)
  cat("Wrote manifest:", manifest_file, "\n")
} else {
  print(manifest_df)
}

cat("Completed.\n")
