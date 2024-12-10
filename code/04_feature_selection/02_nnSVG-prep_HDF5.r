setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(BiocParallel)
	library(nnSVG)
})
set.seed(123)

cat("\nBiocParallel defaults:\n")
MulticoreParam()

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n120_postQC_norm_")
cat("\n\nspe dimensions:\n")
dim(spe)
cat("\ntree:\n")
showtree(logcounts(spe))
cat("\nseed:\n")
seed(logcounts(spe))

spe = filter_genes(spe, filter_genes_ncounts = 3, filter_genes_pcspots = .1, filter_mito=T) #will give 6k genes

spe$dummy_slide = spe$slide
colData(spe)[spe$slide=="V13B23-283","dummy_slide"] = "joint-283-339"
colData(spe)[spe$slide=="V13B23-339","dummy_slide"] = "joint-283-339"

l1 = unique(spe$dummy_slide)
names(l1) = l1
l1 = bplapply(l1, function(x) {
	print(x); Sys.time()
	tmp <- spe[,spe$dummy_slide==x]
	save(tmp, file=paste0("processed-data/04_feature_selection/per-slide_spe/",x,".Rdata"))
})

## Reproducibility information
cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
