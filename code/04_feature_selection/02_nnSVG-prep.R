args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
    stop("Usage: Rscript 02_nnSVG-prep.R <assay_h5>")
}
h5Path <- args[[1L]]
fname = basename(h5Path)
h5_prefix = sub("assays\\.h5$", "", fname)

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(HDF5Array)
    library(DelayedArray)
    library(nnSVG)
})
set.seed(123)
setAutoBlockSize(1e9)


# LOAD NORMALIZED SPE AND FILTER GENES ACROSS ALL SAMPLES #

spe <- loadHDF5SummarizedExperiment(dir=".", prefix=h5_prefix)
stopifnot("sample_id" %in% colnames(colData(spe)), !anyNA(spe$sample_id),
          "logcounts" %in% assayNames(spe))
cat("\nDim spe (before gene filtering):", dim(spe), "\n")
spe <- filter_genes(spe, filter_genes_ncounts=3,
                    filter_genes_pcspots=.1, filter_mito=TRUE)
if (nrow(spe) == 0L) {
    stop("No genes remain after nnSVG preparation filtering.")
}
cat("\nDim spe (filtered genes):", dim(spe), "\n")


# SAVE SELF-CONTAINED PER-SAMPLE HDF5 SPE OBJECTS #

sample_ids <- unique(as.character(spe$sample_id))
dir.create("per-sample_spe", showWarnings=FALSE)
for (sample_id in sample_ids) {
    cat(format(Sys.time()), "Saving sample:", sample_id, "\n")
    tmp <- spe[,spe$sample_id == sample_id]
    saveHDF5SummarizedExperiment(tmp, dir=file.path("per-sample_spe", sample_id),
                               chunkdim=pmin(c(100L,500L), dim(tmp)), verbose=FALSE)
}


# SAVE PER-SAMPLE nnSVG RUN SPECIFICS #

run_specifics <- data.frame(sample_id=sample_ids, nnSVG_seed=123L)
## These two samples produced and irresolvable BiocParallel error that was solved by switching the random seed to 456 in the original analysis.
run_specifics$nnSVG_seed[run_specifics$sample_id %in% c("V13B23-302_B1", "V13B23-333_C1")] <- 456L
run_specifics$n_rows <- NA_integer_
run_specifics$min.genes_per.spot <- NA_integer_
run_specifics$min.spots_per.gene <- NA_integer_
run_specifics$exit_status <- NA_integer_
run_specifics$error_desc <- "Not run"
write.csv(run_specifics, "nnSVG_run_specifics.csv", row.names=FALSE)
