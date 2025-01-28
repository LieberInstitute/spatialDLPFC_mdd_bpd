suppressPackageStartupMessages({
	library(SpatialExperiment)
	#library(scran)
	#library(STexampleData)
	library(spoon)
	#library(BiocParallel)
	library(nnSVG)
	library(DelayedArray)
	library(HDF5Array)
})
set.seed(123)
setAutoBlockSize(1e9)
## for mouse example need to filter to in tissue
#spe <- Visium_mouseCoronal()
#spe <- spe[, colData(spe)$in_tissue == 1]
#spe <- filter_genes(spe)
## for dlpfc example we don't need to filter to in tissue (this example has higher number of spots)
#spe <- spatialLIBD::fetch_data(type="spatialDLPFC_Visium")
#spe = spe[,spe$sample_id=="Br6471_mid"] #4476 in tissue spots
## for dlpfc example less strict filter to get closer to 6k spots
#spe = filter_genes(spe, filter_genes_ncounts = 3, filter_genes_pcspots = .1, filter_mito=T)
## for both mouse and dlpfc reproducible examples recalculate logcounts
#spe <- computeLibraryFactors(spe)
#spe <- logNormCounts(spe)
# for my data
load("processed-data/04_feature_selection/per-slide_spe/V13B23-308.Rdata")
spe = tmp[,tmp$sample_id=="V13B23-308_A1"]
dim(spe)

#nnSVG on DelayedMatrix took forever (2h) so repeat with sparse matrix
cat("\n\nConvert counts and logcounts to dgCMatrix\n")
regular_matrix_counts <- as.matrix(assays(spe)[["counts"]])
sparse_matrix_counts <- as(regular_matrix_counts, "dgCMatrix")
assays(spe)$counts <- sparse_matrix_counts

regular_matrix_lcounts <- as.matrix(assays(spe)[["logcounts"]])
sparse_matrix_lcounts <- as(regular_matrix_lcounts, "dgCMatrix")
assays(spe)$logcounts <- sparse_matrix_lcounts

class(counts(spe))
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

if(class(logcounts(spe))[1]=="DelayedMatrix") weights <- DelayedArray(weights)

cat("\nTime elapsed (min): weighted nnSVG\n")
round(system.time(
	suppressWarnings(output <- weighted_nnSVG(spe, w=weights, n_threads=12))
)[['elapsed']]/60, 2)
rowData(output)

cat("\n\nTime elapsed (min): JT weighted nnSVG\n")
source("code/04_feature_selection/spoon-tutorial_test/weighted-nnSVG_re-write_iterative-cov-matrix.r")
weighted_logcounts <- t(weights)*assays(spe)[['logcounts']]
assay(spe, "weighted_logcounts") <- weighted_logcounts
round(system.time(
        suppressWarnings(output <- weightedJT_nnSVG(spe, X=weights, assay_name="weighted_logcounts", n_threads=12))
)[['elapsed']]/60, 2)
rowData(output)

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
