args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
    stop("Usage: Rscript 03_nnSVG_per-sample.R <sample_spe_dir> <sample_id> <n_threads> <seed>")
}
sample_spe_dir <- args[[1L]]
sample_id <- args[[2L]]
n_threads <- as.integer(args[[3L]])
seed <- as.integer(args[[4L]])
stopifnot(!is.na(n_threads), n_threads > 0L, !is.na(seed))

## Save run specifics in the work directory even when nnSVG fails. 
## Updated throughout the nnSVG process to document progress and failures/errors
run_specifics <- data.frame(
    sample_id=sample_id, nnSVG_seed=seed, n_rows=NA_integer_,
    min.genes_per.spot=NA_integer_, min.spots_per.gene=NA_integer_,
    exit_status=NA_integer_, error_desc="Not completed"
)
write.table(run_specifics, "nnSVG_run_specifics.tsv", sep="\t", row.names=FALSE, quote=FALSE, na="NA")

## Wrap nnSVG per-sample processing and run so that we can specify error behavior within tryCatch (below)
run_nnsvg <- function() {
    suppressPackageStartupMessages({
        library(SpatialExperiment)
        library(HDF5Array)
        library(DelayedArray)
        library(scuttle)
        library(nnSVG)
    })
    set.seed(seed)
    setAutoBlockSize(1e9)


    # LOAD SAMPLE SPE AND FILTER GENES #

    spe <- loadHDF5SummarizedExperiment(dir=sample_spe_dir)
    stopifnot(identical(unique(as.character(spe$sample_id)), sample_id))
    cat(format(Sys.time()), "Sample:", sample_id, "Threads:", n_threads, "Seed:", seed, "\n")
    #cat("\nDim spe (before per-sample gene filtering):", dim(spe), "\n")
    
    ## Following nnSVG guidance, remove genes detected in less than 100 spots
    #### https://bioconductor.org/packages/3.21/bioc/vignettes/nnSVG/inst/doc/nnSVG.html#5_Troubleshooting
    keep_rows = rowSums(counts(spe)>0)>100
    spe = spe[keep_rows,]
    #if (nrow(spe) == 0L) {
    #    stop("No genes are detected in more than 100 spots for sample ", sample_id)
    #}
    spe <- computeLibraryFactors(spe)
    spe <- logNormCounts(spe)
    #cat("\nDim spe (filtered genes):", dim(spe), "\n")
    
    ## Update run specifics with the final gene-filtered, renormalized input passed to nnSVG
    run_specifics$n_rows <<- nrow(spe)
    #### Include min for logcounts (rather than more intuitive counts) because nnSVG take logcounts input and will produce error if not all >0
    run_specifics$min.genes_per.spot <<- min(colSums(logcounts(spe)>0))
    run_specifics$min.spots_per.gene <<- min(rowSums(logcounts(spe)>0))
    write.table(run_specifics, "nnSVG_run_specifics.tsv", sep="\t", row.names=FALSE, quote=FALSE, na="NA")
    #stopifnot(all(colSums(logcounts(spe)) > 0), all(rowSums(logcounts(spe)) > 0))


    # RUN nnSVG AND SAVE RESULTS #

    cat(format(Sys.time()), "Running nnSVG:", sample_id, "\n")
    results <- nnSVG(spe, assay_name="logcounts", n_threads=n_threads)
    write.csv(as.data.frame(rowData(results)), paste0(sample_id, "_nnSVG-results.csv"),
              row.names=TRUE)
    cat(format(Sys.time()), "Saved nnSVG results:", sample_id, "\n")
}

## tryCatch() lets R handle and error with specific code instead of stopping the script,
## allowing us to run nnSVG and catch and document the error for reporting in run specifics
#### These errors should all be resolved but the error process was key to identifying the
#### necessity of using a different random seed for two samples (BiocParallel error)
tryCatch({
    run_nnsvg()
    run_specifics$exit_status <- 0L
    run_specifics$error_desc <- ""
    write.table(run_specifics, "nnSVG_run_specifics.tsv", sep="\t", row.names=FALSE, quote=FALSE, na="NA")
}, error=function(e) {
    run_specifics$exit_status <<- 1L
    # Keep one TSV row; preserve the message while removing tabs and line breaks.
    run_specifics$error_desc <<- gsub("[\t\r\n]+", " ", conditionMessage(e))
    write.table(run_specifics, "nnSVG_run_specifics.tsv", sep="\t", row.names=FALSE, quote=FALSE, na="NA")
    message("nnSVG failed for ", sample_id, ": ", conditionMessage(e))
    quit(save="no", status=1L)
})
