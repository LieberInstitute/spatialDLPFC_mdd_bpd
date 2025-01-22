suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(scran)
	library(STexampleData)
	library(spoon)
	library(BiocParallel)
	library(nnSVG)
})
set.seed(123)

# for mouse example need to filter to in tissue
#spe <- Visium_mouseCoronal()
#spe <- spe[, colData(spe)$in_tissue == 1]
# for dlpfc example we don't need to filter to in tissue (this example has higher number of spots)
spe <- spatialLIBD::fetch_data(type="spatialDLPFC_Visium")
spe = spe[,spe$sample_id=="Br6471_mid"] #4476 in tissue spots
#spe <- filter_genes(spe)
# for dlpfc example less strict filter to get closer to 6k spots
spe = filter_genes(spe, filter_genes_ncounts = 3, filter_genes_pcspots = .1, filter_mito=T)
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)
dim(spe)

# for mouse example further reduce gene set size since so many genes passed nnSVG filter
#cat("\nReduce gene set size\n")
#spe <- spe[sample(rownames(spe), 6061),]
#n10 = c("ENSMUSG00000045691", "ENSMUSG00000029723", "ENSMUSG00000079017",
#        "ENSMUSG00000021816", "ENSMUSG00000025537", "ENSMUSG00000024847",
#        "ENSMUSG00000036620", "ENSMUSG00000026083", "ENSMUSG00000025272",
#        "ENSMUSG00000042729")
#spe <- spe[n10, ]
#dim(spe)

cat("\nAll spots have at least 1 non-zero gene\n")
sum(colSums(logcounts(spe)) > 0)==dim(spe)[2]
cat("All genes have at least 1 non-zero spot\n")
sum(rowSums(logcounts(spe)) > 0)==dim(spe)[1]

cat("\nTime elapsed (min): standard nnSVG\n")
round(system.time(
	nnSVG(spe, assay="logcounts", n_threads=12)
)[['elapsed']]/60, 2)

cat("\nTime elapsed (min): generate weights\n")
round(system.time(
	weights <- generate_weights(input = spe, stabilize = TRUE, n_threads=12)
)[['elapsed']]/60, 2)

cat("\nTime elapsed (min): weighted nnSVG\n")
round(system.time(
	suppressWarnings(output <- weighted_nnSVG(spe, w=weights, n_threads=12))
)[['elapsed']]/60, 2)
rowData(output)

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
