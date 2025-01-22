setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
        library(SpatialExperiment)
        #library(STexampleData)
        library(HDF5Array)
        library(DelayedArray)
        library(scran)
	library(BiocParallel)
        library(nnSVG)
        library(spoon)
})
set.seed(123)

load("processed-data/04_feature_selection/per-slide_spe/V13B23-308.Rdata")
spe <- tmp[,tmp$sample_id=="V13B23-308_A1"]

cat("\n\nConvert logcounts to dgCMatrix",format(Sys.time(), tz="EST"),"...\n")
setAutoBlockSize(1e9)
regular_matrix_counts <- as.matrix(assays(spe)[["counts"]])
sparse_matrix_counts <- as(regular_matrix_counts, "dgCMatrix")
assays(spe)$counts <- sparse_matrix_counts
cat(format(Sys.time(), tz="EST"))

cat("\n\nFilter genes and re-log norm counts",format(Sys.time(), tz="EST"),"...\n")
spe <- filter_genes(spe)
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)


#weights reproducibility test
weightList <- list()
## from matrix
cat("\n\nGenerate weights from matrix",format(Sys.time(), tz="EST"),"...\n")
set.seed(123)
weightList$from_mtx <- generate_weights(input = logcounts(spe), 
	spatial_coords=spatialCoords(spe), stabilize = TRUE, n_threads=12)
## from spe
cat("\n\nGenerate weights from spe",format(Sys.time(), tz="EST"),"...\n")
set.seed(123)
weightList$from_spe <- generate_weights(input = spe, stabilize = TRUE, n_threads=12)
## checkpoint save
cat("\n\nCheckpoint save weightList",format(Sys.time(), tz="EST"),"...\n")
saveRDS(weightList, "processed-data/04_feature_selection/spoon-tutorial_test/weights_reproducibility_V13B23-308_A1.rda")
## with rngseed
cat("\n\nGenerate weights from matrix with RNG seed",format(Sys.time(), tz="EST"),"...\n")
cat("(Format used in 'code/04_feature_selection/test_spoon-nnSVG_per-sample.r')\n")
set.seed(123)
weightList$RNG_seed4 <- generate_weights(input = logcounts(spe),
        spatial_coords=spatialCoords(spe), stabilize = TRUE, #n_threads=12,
	BPPARAM = MulticoreParam(workers=12, RNGseed = 4))

cat("\n\nSave completed weightList",format(Sys.time(), tz="EST"),"...\n")
saveRDS(weightList, "processed-data/04_feature_selection/spoon-tutorial_test/weights_reproducibility_V13B23-308_A1.rda")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
