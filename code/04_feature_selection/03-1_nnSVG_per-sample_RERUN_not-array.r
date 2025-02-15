suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(HDF5Array)
  library(DelayedArray)
  library(scran)
  library(nnSVG)
})
set.seed(123)

x = readLines("processed-data/04_feature_selection/per-sample_spe_RERUN_list.txt")

cat(x[1])
load(paste0("processed-data/04_feature_selection/per-sample_spe/",x[1]))
regular_matrix_counts <- as.matrix(assays(tmp)[["counts"]])
sparse_matrix_counts <- as(regular_matrix_counts, "dgCMatrix")
assays(tmp)$counts <- sparse_matrix_counts

keep_rows = rowSums(counts(tmp))!=0
table(keep_rows)
tmp = tmp[keep_rows,]
tmp <- computeLibraryFactors(tmp)
tmp <- logNormCounts(tmp)

dim(tmp)
name1 = substr(x[1], start=0, stop=13)
cat("\n\nStart nnSVG",name1,format(Sys.time(), tz="EST"),"\n\n")
results <- nnSVG(tmp, assay_name="logcounts", n_threads=12, verbose=F)
cat("\n\nFinished",name1,format(Sys.time(),tz="EST"),"\n\n")
write.csv(rowData(results),"processed-data/04_feature_selection/per-sample_svgs/",name1,"_nnSVG-results.csv",row.names=T)


cat("\n\n",x[2])
load(paste0("processed-data/04_feature_selection/per-sample_spe/",x[2]))
regular_matrix_counts <- as.matrix(assays(tmp)[["counts"]])
sparse_matrix_counts <- as(regular_matrix_counts, "dgCMatrix")
assays(tmp)$counts <- sparse_matrix_counts

keep_rows = rowSums(counts(tmp))!=0
table(keep_rows)
tmp = tmp[keep_rows,]
tmp <- computeLibraryFactors(tmp)
tmp <- logNormCounts(tmp)

dim(tmp)
name2 = substr(x[2], start=0, stop=13)
cat("\n\nStart nnSVG",name2,format(Sys.time(), tz="EST"),"\n\n")
results <- nnSVG(tmp, assay_name="logcounts", n_threads=12, verbose=F)
cat("\n\nFinished",name2,format(Sys.time(),tz="EST"),"\n\n")
write.csv(rowData(results),"processed-data/04_feature_selection/per-sample_svgs/",name2,"_nnSVG-results.csv",row.names=T)


cat("\n\n",x[3])
load(paste0("processed-data/04_feature_selection/per-sample_spe/",x[3]))
regular_matrix_counts <- as.matrix(assays(tmp)[["counts"]])
sparse_matrix_counts <- as(regular_matrix_counts, "dgCMatrix")
assays(tmp)$counts <- sparse_matrix_counts

keep_rows = rowSums(counts(tmp))!=0
table(keep_rows)
tmp = tmp[keep_rows,]
tmp <- computeLibraryFactors(tmp)
tmp <- logNormCounts(tmp)

dim(tmp)
name3 =	substr(x[3], start=0, stop=13)
cat("\n\nStart nnSVG",name3,format(Sys.time(), tz="EST"),"\n\n")
results <- nnSVG(tmp, assay_name="logcounts", n_threads=12, verbose=F)
cat("\n\nFinished",name3,format(Sys.time(),tz="EST"),"\n\n")
write.csv(rowData(results),"processed-data/04_feature_selection/per-sample_svgs/",name3,"_nnSVG-results.csv",row.names=T)


cat("\n\n",x[4])
load(paste0("processed-data/04_feature_selection/per-sample_spe/",x[4]))
regular_matrix_counts <- as.matrix(assays(tmp)[["counts"]])
sparse_matrix_counts <- as(regular_matrix_counts, "dgCMatrix")
assays(tmp)$counts <- sparse_matrix_counts

keep_rows = rowSums(counts(tmp))!=0
table(keep_rows)
tmp = tmp[keep_rows,]
tmp <- computeLibraryFactors(tmp)
tmp <- logNormCounts(tmp)

dim(tmp)
name4 =	substr(x[4], start=0, stop=13)
cat("\n\nStart nnSVG",name2,format(Sys.time(), tz="EST"),"\n\n")
results <- nnSVG(tmp, assay_name="logcounts", n_threads=12, verbose=F)
cat("\n\nFinished",name4,format(Sys.time(),tz="EST"),"\n\n")
write.csv(rowData(results),"processed-data/04_feature_selection/per-sample_svgs/",name4,"_nnSVG-results.csv",row.names=T)


## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
