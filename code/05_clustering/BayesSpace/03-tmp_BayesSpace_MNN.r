setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(HDF5Array)
  library(DelayedArray)
  library(BayesSpace)
})
set.seed(123)

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n120_postQC_norm_")
if(sum(c("row","col") %in% colnames(colData(spe)))!=2) {
	spe$row <- spe$array_row
	spe$col <- spe$array_col
}

## -- 2026-08-13 update --
## Note that the spatial coordinates were not offset and should have been.
## See http://edward130603.github.io/BayesSpace/articles/joint_clustering.html#clustering-1
## for more details.

load("processed-data/05_clustering/BayesSpace/spe_fastMNN_n1079-k10-ordered-merge.Rdata")
stopifnot(identical(rownames(mnn), rownames(colData(spe))))
reducedDim(spe, "MNN_1079") <- mnn$corrected
#temporary code while PCs finish running
#load("processed-data/05_clustering/BayesSpace/spe_PCA_MNN_empty-assays.Rdata")
#reducedDim(spe, "PCA_1079") <- reducedDim(dummy_spe, "PCA_1079")

spe

cat("\nspatialPreprocess:",format(Sys.time()),"\n")
spe <- spatialPreprocess(spe, platform="Visium", skip.PCA=TRUE, log.normalize=FALSE)
#spe <- spatialPreprocess(spe, platform="Visium", skip.PCA=FALSE, log.normalize=FALSE,
#	n.PCs=15, n.HVGs=2000)

#mimic the labeling of top HVGs performed by default spatial preprocess in case it is used downstream
#rowData(spe)$is.HVG = rownames(spe) %in% rownames(attr(reducedDim(dummy_spe, "PCA_1079"),"rotation"))

cat("\n\nspatialCluster:",format(Sys.time()),"\n")
(params = "BayesSpace_MNN.1079.d13_q9_iter.1k_init.PRECAST")
spe <- spatialCluster(spe, q=9, use.dimred="MNN_1079", d=13,
	init=spe$precast_k9_1079,
	#init.method="kmeans", #initial clustering method options = mclust or kmeans
	model="t", #Using t-distributed errors in likelihood/bayesian models is a different way to obtain robust methods, as the t-distribution has heavier tails than the normal.
	nrep=1000, #default nrep is 50k, vignette recommends min 10k so will start there
	burn.in=100
)

#isolate result columns
tmp = colData(spe)[,c("row","col","spot.idx","spot.neighbors","cluster.init","spatial.cluster")]
#rename spatial.cluster to be more specific
colnames(tmp)[6] = params
#save results as csv
write.csv(tmp, paste0("processed-data/05_clustering/BayesSpace/colData_",gsub("\\.","-",params),".csv"), row.names=T)


#colnames(colData(spe))[grep("spatial", colnames(colData(spe)))] = "BayesSpace_PCA.1079.d20_q9"
#table(colData(spe)[,c("precast_k9_1079","BayesSpace_PCA.1079.d20_q9")])

#quickResaveHDF5SummarizedExperiment(spe)
#cat("\nSaved with quickResave to spe HDF5 dir: processed-data/04_feature_selection/spe_n120_postQC_norm_\n")
#dummy_spe = spe
#assays(dummy_spe, "counts") <- NULL
#assays(dummy_spe, "logcounts") <- NULL
#cat("\n\nSave results:",format(Sys.time()),"\n")
#save(dummy_spe, file=paste0("processed-data/05_clustering/BayesSpace/spe_",gsub("\\.","-",params),"_empty-assays.Rdata"))
#cat("Saved to:",paste0("processed-data/05_clustering/spe_",gsub("\\.","-",params),"_empty-assays.Rdata"))

## Reproducibility information
cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
