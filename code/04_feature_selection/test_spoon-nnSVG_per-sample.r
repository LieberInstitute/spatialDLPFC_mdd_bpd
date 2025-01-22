args = commandArgs(TRUE)
print(args[[1]])
setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(scran)
	library(spoon)
	library(nnSVG)
	library(HDF5Array)
	library(DelayedArray)
	library(BiocParallel)
})
### rather than parallelizing the nnSVG() function over each gene, using the gene-specific weights matrix column as the covariate, (Kinnary's weighted_nnSVG)
### instead modify the BRISC parallelization within the nnSVG() function so that only the column of the weight matrix corresponding to the relevant gene is used as a covariate
source("code/04_feature_selection/nnSVG_re-write_iterative-cov-matrix.r")
set.seed(123)

load("processed-data/04_feature_selection/per-slide_spe/V13B23-308.Rdata")
dim(tmp)
tmp = tmp[,tmp$sample_id==args[[1]]]
dim(tmp)

cat("\n\nFilter out zero genes (if present) and generate logcounts on subset of genes...\n")
keep_rows = rowSums(counts(tmp))!=0
tmp = tmp[keep_rows,]
tmp <- computeLibraryFactors(tmp)
tmp <- logNormCounts(tmp)
dim(logcounts(tmp))

cat("\n\nConvert logcounts to dgCMatrix",format(Sys.time(), tz="EST"),"...\n")
setAutoBlockSize(1e9)
regular_matrix_logcounts <- as.matrix(assays(tmp)[["logcounts"]])
sparse_matrix_logcounts <- as(regular_matrix_logcounts, "dgCMatrix")
assays(tmp)$logcounts <- sparse_matrix_logcounts
cat(format(Sys.time(), tz="EST"))
colnames(logcounts(tmp))[1:4]
rownames(logcounts(tmp))[1:4]

#if weights exist then load them, else calculate
#if(paste0("test_spoon-nnSVG_per-sample_weights-",args[[1]],".csv") %in% list.files("processed-data/04_feature_selection")) {
#	cat("\n\nLoading saved weights...")
#	weights = read.csv(paste0("processed-data/04_feature_selection/test_spoon-nnSVG_per-sample_weights-",args[[1]],".csv"), row.names=1)
#}
#else {
	cat("\n\nGenerate weights", format(Sys.time(), tz="EST"),"...\n")
	weights = generate_weights(input=logcounts(tmp), spatial_coords=spatialCoords(tmp), stabilize=T,
		BPPARAM = MulticoreParam(workers = 12, RNGseed = 4))
	cat(format(Sys.time(), tz="EST"))
	write.csv(weights, paste0("processed-data/04_feature_selection/test_spoon-nnSVG_per-sample_weights-",args[[1]],".csv"), row.names=T)
	cat("\nWeights saved to:", paste0("processed-data/04_feature_selection/test_spoon-nnSVG_per-sample_weights-",args[[1]],".csv"))
#}

if(dim(weights)[1]==dim(logcounts(tmp))[2]) {
	rownames(weights) = colnames(logcounts(tmp))
}
if(dim(weights)[2]==dim(logcounts(tmp))[1]) {
	colnames(weights) = rownames(logcounts(tmp))
}

# calculate weighted logcount matrix
### from here: https://github.com/kinnaryshah/spoon/blob/main/R/weighted_nnSVG.R#L109
weighted_logcounts <- t(weights)*assays(tmp)[['logcounts']]
weighted_mean <- Matrix::rowMeans(weighted_logcounts)
assay(tmp, "weighted_logcounts") <- weighted_logcounts

#weighted nnSVG
cat("\n\nWeighted nnSVG",format(Sys.time(), tz="EST"),"...")
output = weightedJT_nnSVG(input=tmp, assay="weighted_logcounts", X=weights,
	n_threads=12)
cat("\n")
cat(format(Sys.time(), tz="EST"))
write.csv(rowData(output), paste0("processed-data/04_feature_selection/test_spoon-nnSVG_per-sample_results-",args[[1]],".csv"), row.names=T)
cat("\nWeighted nnSVG results saved to:", paste0("processed-data/04_feature_selection/test_spoon-nnSVG_per-sample_results-",args[[1]],".csv"))



## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
