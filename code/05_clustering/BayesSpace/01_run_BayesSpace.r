setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(HDF5Array)
  library(DelayedArray)
  library(BayesSpace)
})
set.seed(123)

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n120_postQC_norm_")
spe$row <- spe$array_row
spe$col <- spe$array_col

load("processed-data/04_feature_selection/spe_PCA_empty-assays.Rdata")
reducedDim(spe, "PCA_1054") <- reducedDim(dummy_spe, "PCA_1054")

spe

cat("\nspatialPreprocess:",format(Sys.time()),"\n")
spe <- spatialPreprocess(spe, platform="Visium", skip.PCA=TRUE, log.normalize=FALSE)

cat("\n\nspatialCluster:",format(Sys.time()),"\n")
spe <- spatialCluster(spe, q=7, use.dimred="PCA_1054", d=20, #start with using 20 MNN dims
	init.method="kmeans", #initial clustering method options = mclust or kmeans
	model="t", #Using t-distributed errors in likelihood/bayesian models is a different way to obtain robust methods, as the t-distribution has heavier tails than the normal.
	nrep=10000 #default nrep is 50k, vignette recommends min 10k so will start there
)
table(colData(spe)[,"spatial.cluster"])

#rename results more descriptive
colnames(colData(spe))[grep("spatial", colnames(colData(spe)))] = "BayesSpace_PCA.1054_q7"
table(colData(spe)[,"BayesSpace_PCA.1054_k7"])

dummy_spe = spe
assays(dummy_spe, "counts") <- NULL
assays(dummy_spe, "logcounts") <- NULL
cat("\n\nSave results:",format(Sys.time()),"\n")
save(dummy_spe, file="processed-data/05_clustering/spe_bayes_PCA-1054_k7-kmeans-10k_empty-assays.Rdata")
cat("Saved to: processed-data/05_clustering/spe_bayes_PCA-1054_k7-kmeans-10k_empty-assays.Rdata")

## Reproducibility information
cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
