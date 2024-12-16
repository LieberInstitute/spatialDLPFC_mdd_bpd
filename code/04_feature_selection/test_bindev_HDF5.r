setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
        library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
        library(scry)
        library(dplyr)
})
set.seed(123)

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n120_postQC_norm_")
#read in any nnSVG result to get the feature list after filtering
tmp = read.csv("processed-data/04_feature_selection/per-slide_svgs/V13B23-301_nnSVG-results.csv", row.names=1)
spe <- spe[rownames(tmp),]
dim(spe)
cat("\nRealize subset spe...\n")
format(Sys.time(), tz="EST")
rlz=writeHDF5Array(counts(spe)) 
getHDF5DumpDir()
format(Sys.time(), tz="EST")
dim(rlz)

cat("\nRunning default model...\n")
format(Sys.time(), tz="EST")
default <- devianceFeatureSelection(rlz, fam="binomial", batch=NULL)

df = cbind.data.frame(rowData(spe)[names(default),c("gene_id","gene_name")],"dev"=default,
        "rank"=(length(default)+1)-rank(default))

cat("\nRunning batch model...\n")
format(Sys.time(), tz="EST")
batch <- devianceFeatureSelection(rlz, fam="binomial", batch=as.factor(spe$slide))

df = left_join(df, cbind.data.frame(rowData(spe)[rownames(rlz),c("gene_id","gene_name")], "dev"=batch, "rank"=(length(batch)+1)-rank(batch)),
	by=c("gene_id","gene_name"), suffix=c("_default","_batch")) 

cat("\nSave output...\n")
format(Sys.time(), tz="EST")
write.csv(df, "processed-data/04_feature_selection/test_bindev_batch-slide.csv",row.names=FALSE)
cat("\nSaved to: processed-data/04_feature_selection/test_bindev_batch-slide.csv\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
