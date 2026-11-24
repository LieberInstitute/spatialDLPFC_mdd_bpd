args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
    stop("Usage: Rscript 01_apply-QC_norm.R <assay_h5> <qc_coldata>")
}
h5Path <- args[[1L]]
fname = basename(h5Path)
h5_prefix = sub("assays\\.h5$", "", fname)
output_prefix <- paste0(h5_prefix, "postQC_norm_")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(HDF5Array)
    library(DelayedArray)
    library(scuttle)
})
set.seed(123)
setAutoBlockSize(1e9)


# LOAD SPE AND QC ANNOTATIONS #

qc.df = read.csv(args[[2L]], row.names=1)
spe <- loadHDF5SummarizedExperiment(dir=".", prefix=h5_prefix)
stopifnot(identical(rownames(qc.df), colnames(spe)),
          "remove_spots" %in% colnames(qc.df),
          is.logical(qc.df$remove_spots), !anyNA(qc.df$remove_spots))
cat("\nDim spe (before QC):", dim(spe), "\n")


# APPLY QC FILTERS #

spe$remove_spots = qc.df$remove_spots
spe = spe[,!spe$remove_spots]
cat("\nDim spe (QC filtered):", dim(spe), "\n")

# Stage 1 already excludes the extra V13F27-338_C1 sample from the n119 dataset.
spots_more0 = rowSums(counts(spe)>0)
rowData(spe)$n_spots = spots_more0
spe = spe[spots_more0>5,]
cat("\nDim spe (genes detected in <=5 spots removed):", dim(spe), "\n")


# NORMALIZATION #

cat(format(Sys.time()), "Compute library factors and normalize counts\n")
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)


# SAVE FILTERED AND NORMALIZED SPE #

start.time = Sys.time()
cat("\nStart HDF5SummarizedExperiment save:", format(start.time), "\n")
saveHDF5SummarizedExperiment(spe, dir=".", prefix=output_prefix,
                           chunkdim=pmin(c(100L,500L), dim(spe)), verbose=FALSE)
cat("\nTime elapsed (saveHDF5):",
    round(difftime(Sys.time(), start.time, units="hours"), 2), "hours\n")
