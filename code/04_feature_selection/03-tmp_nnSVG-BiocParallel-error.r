suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(HDF5Array)
  library(DelayedArray)
  library(scran)
  library(nnSVG)
})
set.seed(123)

load("processed-data/04_feature_selection/per-sample_spe/V13B23-279_A1.Rdata")
a1 = tmp
regular_matrix_counts <- as.matrix(assays(a1)[["counts"]])
sparse_matrix_counts <- as(regular_matrix_counts, "dgCMatrix")
assays(a1)$counts <- sparse_matrix_counts

keep_rows = rowSums(counts(a1))!=0
table(keep_rows)
a1 = a1[keep_rows,]
a1 <- computeLibraryFactors(a1)
a1 <- logNormCounts(a1)

dim(a1)

load("processed-data/04_feature_selection/per-sample_spe/V13B23-279_C1.Rdata")
c1 = tmp
regular_matrix_counts <- as.matrix(assays(c1)[["counts"]])
sparse_matrix_counts <- as(regular_matrix_counts, "dgCMatrix")
assays(c1)$counts <- sparse_matrix_counts

keep_rows = rowSums(counts(c1))!=0
table(keep_rows)
c1 = c1[keep_rows,]
c1 <- computeLibraryFactors(c1)
c1 <- logNormCounts(c1)

dim(c1)

#run in interactive session
##compare spots per gene
#spots_per_gene = rowSums(logcounts(a1)>0)
#spots_per_gene2 = rowSums(logcounts(c1)>0)

#plot(ecdf(spots_per_gene), col="red", main="Failed in red; Comparison in black", xlab="# of spots with non-0 logcounts per gene")
#lines(ecdf(spots_per_gene2))

cat("\n\nStart nnSVG V13B23-279_C1 (completed)...",format(Sys.time(), tz="EST"),"\n\n")
results_c1 <- nnSVG(c1, assay_name="logcounts", n_threads=12, verbose=T)
cat("\n\nFinished",format(Sys.time(),tz="EST"),"\n\n")

cat("\n\nStart nnSVG V13B23-279_A1 (BiocParallel error)...",format(Sys.time(), tz="EST"),"\n\n")
results_a1 <- nnSVG(a1, assay_name="logcounts", n_threads=12, verbose=T)
cat("\n\nFinished",format(Sys.time(),tz="EST"),"\n")
write.csv(rowData(results_a1), "processed-data/04_feature_selection/per-sample_svgs/V13B23-279_A1_nnSVG-results.csv", row.names=T)
cat("saved results to: processed-data/04_feature_selection/per-sample_svgs/V13B23-279_A1_nnSVG-results.csv")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
