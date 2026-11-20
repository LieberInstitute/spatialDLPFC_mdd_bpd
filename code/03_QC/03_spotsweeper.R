args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(3L, 4L)) {
    stop("Usage: Rscript 03_spotsweeper.R <spe_dir> <n_samples> <qc_coldata> [comma-separated_ignore_edge_samples]")
}
spe_dir <- args[[1L]]
n_samples <- as.integer(args[[2L]])
spe_prefix <- paste0("spe_n", n_samples, "_")
ignore_edge_samples <- if (length(args) == 4L && nzchar(args[[4L]])) {
    strsplit(args[[4L]], ",", fixed=TRUE)[[1L]]
} else {
    character()
}

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(HDF5Array)
    library(SpotSweeper)
})


# LOAD SPE AND QC ANNOTATIONS #

spe <- loadHDF5SummarizedExperiment(dir=spe_dir, prefix=spe_prefix)
stopifnot(length(unique(spe$sample_id)) == n_samples)

qc.df = read.csv(args[[3L]], row.names=1)
stopifnot(identical(colnames(spe), rownames(qc.df)))
required_columns <- c("sample_id", "sum_umi", "edge_outlier", "outlier_area")
stopifnot(all(required_columns %in% colnames(qc.df)))
if (!all(ignore_edge_samples %in% unique(qc.df$sample_id))) {
    stop("Unknown ignore_edge_samples: ",
         paste(setdiff(ignore_edge_samples, unique(qc.df$sample_id)), collapse=", "))
}

## If provided, ignore edge label for specific samples
if (length(ignore_edge_samples) > 0L) {
    message("Ignoring edge labels for samples: ",
            paste(unique(ignore_edge_samples), collapse=", "))
} else {
    message("Keeping default edge labels for all samples")
}
qc.df$true_edges = ifelse(qc.df$sample_id %in% c(ignore_edge_samples), FALSE, qc.df$edge_outlier)

qc.df$lowUMI = qc.df$sum_umi <= 100
qc.df$remove_spots = qc.df$true_edges | qc.df$outlier_area | qc.df$lowUMI

## Transfer QC annotations to SPE and subset for spotsweeper
spe$remove_spots = qc.df$remove_spots
spe = spe[,!spe$remove_spots & spe$in_tissue]



# RUN SPOTSWEEPER #

## Calculate local outliers: UMI
spe <- localOutliers(spe, metric = "sum_umi", direction = "lower", log = TRUE)
colnames(colData(spe))[ncol(colData(spe))-1] = "umi_local.outlier"

## Calculate local outliers: detected genes
spe <- localOutliers(spe, metric = "sum_gene", direction = "lower", log = TRUE)
colnames(colData(spe))[ncol(colData(spe))-1] = "genes_local.outlier"

## Calculate local outliers: mito fractions
spe <- localOutliers(spe, metric = "expr_chrM_ratio", direction = "higher", log = FALSE)
colnames(colData(spe))[ncol(colData(spe))-1] = "chrM.ratio_local.outlier"

# UPDATE AND SAVE QC CSV #
true.spots = colnames(spe)

qc.df$umi_local.outlier = NA
qc.df[true.spots, "umi_local.outlier"] = colData(spe)$umi_local.outlier
qc.df$genes_local.outlier = NA
qc.df[true.spots, "genes_local.outlier"] = colData(spe)$genes_local.outlier
qc.df$chrM.ratio_local.outlier = NA
qc.df[true.spots, "chrM.ratio_local.outlier"] =	colData(spe)$chrM.ratio_local.outlier

qc.df$spotsweeper_outlier = qc.df$umi_local.outlier | qc.df$genes_local.outlier | qc.df$chrM.ratio_local.outlier

qc.df$remove_spots = qc.df$remove_spots | qc.df$spotsweeper_outlier
c1 = grep("remove_spots", colnames(qc.df))
qc.df = qc.df[,c(setdiff(1:ncol(qc.df), c1), c1)]

#save coldata
write.csv(qc.df, "colData_QC.csv", row.names=T)
