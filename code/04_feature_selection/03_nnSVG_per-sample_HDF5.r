args = commandArgs(TRUE)
print(args[[1]])
setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(HDF5Array)
  library(DelayedArray)
  library(scran)
  library(nnSVG)
})
set.seed(123)

load(file=paste0("processed-data/04_feature_selection/per-slide_spe/",args[[1]]))
dim(tmp)

cat("\n\nConvert to dgCMatrix...")
format(Sys.time(), tz="EST")
regular_matrix_counts <- as.matrix(assays(tmp)[["counts"]])
sparse_matrix_counts <- as(regular_matrix_counts, "dgCMatrix")
assays(tmp)$counts <- sparse_matrix_counts

cat("\nSubsetting to individual samples...\n")
format(Sys.time(), tz="EST")
l1 = unique(tmp$sample_id)
names(l1) = l1
spe.list = lapply(l1, function(x) tmp[,tmp$sample_id==x])

cat("\nLooping nnSVG...\n")
for(i in seq_along(spe.list)) {
	cat(names(spe.list)[i])
	cat("\nFilter out zero genes (if present) and generate logcounts on subset of genes...\n")
	keep_rows = rowSums(counts(spe.list[[i]]))!=0
	tmp_sub = spe.list[[i]][keep_rows,]
	tmp_sub <- computeLibraryFactors(tmp_sub)
	tmp_sub <- logNormCounts(tmp_sub)
	cat(dim(logcounts(tmp_sub)))
	cat("\nnnSVG start...",format(Sys.time(),tz="EST"),"\n")
	set.seed(123)
	results <- nnSVG(tmp_sub, n_threads=12)
	svg = rowData(results)
	cat(names(spe.list)[i],"save output",format(Sys.time(),tz="EST"),"\n\n")
	write.csv(svg, paste0("processed-data/04_feature_selection/per-sample_svgs/", names(spe.list)[i], "_nnSVG-results.csv"), row.names=T)
}



## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
