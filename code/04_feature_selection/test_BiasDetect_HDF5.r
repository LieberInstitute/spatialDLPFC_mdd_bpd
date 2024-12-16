setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
        library(SpatialExperiment)
        library(HDF5Array)
        library(DelayedArray)
        library(BiasDetect)
})
set.seed(123)

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n120_postQC_norm_")
dim(spe)
#read in any nnSVG result to get the feature list after filtering
tmp = read.csv("processed-data/04_feature_selection/per-slide_svgs/V13B23-301_nnSVG-results.csv", row.names=1)
#cat("n genes to test:", length(rownames(tmp)))
spe = spe[rownames(tmp),]
dim(spe)

spe$dummy_slide = as.character(spe$slide)
colData(spe)[spe$slide=="V13B23-283","dummy_slide"] = "joint-283-339"
colData(spe)[spe$slide=="V13B23-339","dummy_slide"] = "joint-283-339"


cat("\n\nTry on subset spe...\n")
format(Sys.time(), tz="EST")
batch_df = featureSelect(spe, batch="dummy_slide", VGs=rownames(tmp))
format(Sys.time(), tz="EST")
write.csv(batch_df, "processed-data/04_feature_selection/test_BiasDetect-featureSelect_dummy-slide.csv", row.names=F)
cat("\nSaved batch DF to: processed-data/04_feature_selection/test_BiasDetect-featureSelect_dummy-slide.csv\n")

outList = list()
outList$plots = biasDetect(batch_df, nSD_dev=5, visual=T)
outList$genes = biasDetect(batch_df, nSD_dev=5, visual=F)
saveRDS(outList, "processed-data/04_feature_selection/test_BiasDetect-biasDetect.rda")
cat("\nSaved plots and biased gene list to: processed-data/04_feature_selection/test_BiasDetect-biasDetect_dummy-slide.rda\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
