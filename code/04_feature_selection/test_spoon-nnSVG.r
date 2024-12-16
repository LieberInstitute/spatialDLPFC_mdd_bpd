setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(spoon)
	library(nnSVG)
	library(HDF5Array)
	library(DelayedArray)
	library(BiocParallel)
})
set.seed(123)

(i="joint-283-339")
load(paste0("processed-data/04_feature_selection/per-slide_spe/",i,".Rdata"))
dim(tmp)

cat("\n\nRealize logcounts matrix...")
format(Sys.time(), tz="EST")
rlz = writeHDF5Array(logcounts(tmp))
dim(rlz)
showtree(rlz)
seed(rlz)
format(Sys.time(), tz="EST")

cat("\n\nGenerate weights...")
format(Sys.time(), tz="EST")
weights = generate_weights(input=rlz, spatial_coords=spatialCoords(tmp), stabilize=T,
	BPPARAM = MulticoreParam(workers = 12, RNGseed = 4))
format(Sys.time(), tz="EST")
write.csv(weights, "processed-data/04_feature_selection/test_spoon-nnSVG_weights.csv", row.names=T)
cat("\n\nWeights saved to: processed-data/04_feature_selection/test_spoon-nnSVG_weights.csv")

cat("\n\nWeighted nnSVG...")
format(Sys.time(), tz="EST")
output = weighted_nnSVG(input=rlz, w=weights, spatial_coords=spatialCoords(tmp), 
	BPPARAM = MulticoreParam(workers=12, RNGseed = 5))
format(Sys.time(), tz="EST")

write.csv(output, "processed-data/04_feature_selection/test_spoon-nnSVG_results.csv", row.names=T)
cat("\n\nWeighted nnSVG results saved to: processed-data/04_feature_selection/test_spoon-nnSVG_results.csv") 

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
