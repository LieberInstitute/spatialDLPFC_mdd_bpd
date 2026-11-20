args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2L, 3L)) {
    stop("Usage: Rscript 03_spotsweeper.R <assay_h5> <qc_df> [comma-separated_ignore_edge_samples]")
}
h5Path <- args[[1L]]
fname = basename(h5Path)
h5_prefix = sub("assays\\.h5$", "", fname)

ignore_edge_samples <- if (length(args) == 3L && nzchar(args[[3L]])) {
    strsplit(args[[3L]], ",", fixed=TRUE)[[1L]]
} else {
    character()
}

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(HDF5Array)
    library(SpotSweeper)
})


# LOAD SPE AND QC ANNOTATIONS #

spe <- loadHDF5SummarizedExperiment(dir=".", prefix=h5_prefix)

qc.df = read.csv(args[[2L]], row.names=1)
stopifnot(identical(colnames(spe), rownames(qc.df)))

## If provided, ignore edge label for specific samples
if (!all(ignore_edge_samples %in% unique(qc.df$sample_id))) {
    warning("Provided ignore_edge_samples are not in the current subset: ",
            paste(setdiff(ignore_edge_samples, unique(qc.df$sample_id)), collapse=", "),
            ". Continuing with samples present in the subset.", call.=FALSE)
}
ignore_edge_samples <- intersect(ignore_edge_samples, unique(qc.df$sample_id))

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

qc.df$umi_local.outlier = FALSE
qc.df[true.spots, "umi_local.outlier"] = colData(spe)$umi_local.outlier
qc.df$genes_local.outlier = FALSE
qc.df[true.spots, "genes_local.outlier"] = colData(spe)$genes_local.outlier
qc.df$chrM.ratio_local.outlier = FALSE
qc.df[true.spots, "chrM.ratio_local.outlier"] =	colData(spe)$chrM.ratio_local.outlier

qc.df$spotsweeper_outlier = qc.df$umi_local.outlier | qc.df$genes_local.outlier | qc.df$chrM.ratio_local.outlier

qc.df$remove_spots = qc.df$remove_spots | qc.df$spotsweeper_outlier | !qc.df$in_tissue
c1 = grep("remove_spots", colnames(qc.df))
qc.df = qc.df[,c(setdiff(1:ncol(qc.df), c1), c1)]

#save coldata
write.csv(qc.df, "colData_QC.csv", row.names=T)
