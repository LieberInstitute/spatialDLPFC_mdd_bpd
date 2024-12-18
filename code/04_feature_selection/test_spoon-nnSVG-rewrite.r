setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	#library(spoon)
	library(nnSVG)
	library(HDF5Array)
	library(DelayedArray)
	library(BiocParallel)
})
set.seed(123)

(i="joint-283-339")
load(paste0("processed-data/04_feature_selection/per-slide_spe/",i,".Rdata"))
dim(tmp)

#cat("\n\nRealize logcounts matrix...")
#format(Sys.time(), tz="EST")
#rlz = writeHDF5Array(logcounts(tmp))
#dim(rlz)
#showtree(rlz)
#seed(rlz)
#format(Sys.time(), tz="EST")

cat("\n\nConvert to dgCMatrix...")
format(Sys.time(), tz="EST")
regular_matrix_logcounts <- as.matrix(assays(tmp)[["logcounts"]])
sparse_matrix_logcounts <- as(regular_matrix_logcounts, "dgCMatrix")
assays(tmp)$logcounts <- sparse_matrix_logcounts
format(Sys.time(), tz="EST")
colnames(logcounts(tmp))[1:4]
rownames(logcounts(tmp))[1:4]
#cat("\n\nGenerate weights...")
#format(Sys.time(), tz="EST")
#weights = generate_weights(input=rlz, spatial_coords=spatialCoords(tmp), stabilize=T,
#	BPPARAM = MulticoreParam(workers = 12, RNGseed = 4))
#format(Sys.time(), tz="EST")
#write.csv(weights, "processed-data/04_feature_selection/test_spoon-nnSVG_weights.csv", row.names=T)
#cat("\n\nWeights saved to: processed-data/04_feature_selection/test_spoon-nnSVG_weights.csv")

cat("\n\nLoad saved weights...")
format(Sys.time(), tz="EST")
weights = read.csv("processed-data/04_feature_selection/test_spoon-nnSVG_weights.csv", row.names=1)
format(Sys.time(), tz="EST")

if(dim(weights)[1]==dim(logcounts(tmp))[2]) {
	rownames(weights) = colnames(logcounts(tmp))
}
if(dim(weights)[2]==dim(logcounts(tmp))[1]) {
	colnames(weights) = rownames(logcounts(tmp))
}


cat("\n\nWeighted nnSVG...")
format(Sys.time(), tz="EST")
# calculate weighted logcount matrix
### from here: https://github.com/kinnaryshah/spoon/blob/main/R/weighted_nnSVG.R#L109
weighted_logcounts <- t(weights)*assays(tmp)[['logcounts']]
weighted_mean <- Matrix::rowMeans(weighted_logcounts)
assay(tmp, "weighted_logcounts") <- weighted_logcounts

# now try my re-written nnSVG code
### rather than parallelizing the nnSVG() function over each gene, using the gene-specific weights matrix column as the covariate, (Kinnary's weighted_nnSVG)
### instead modify the BRISC parallelization within the nnSVG() function so that only the column of the weight matrix corresponding to the relevant gene is used as a covariate
source("code/04_feature_selection/nnSVG_re-write_iterative-cov-matrix.r")
output = weighted_nnSVG(input=tmp, assay="weighted_logcounts", X=weights, 
	n_threads=12)
format(Sys.time(), tz="EST")

write.csv(rowData(output), "processed-data/04_feature_selection/test_spoon-nnSVG-rewrite_results.csv", row.names=T)
cat("\n\nWeighted nnSVG results saved to: processed-data/04_feature_selection/test_spoon-nnSVG-rewrite_results.csv") 

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
