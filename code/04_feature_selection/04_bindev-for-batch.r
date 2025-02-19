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
spe$dummy_slide = ifelse(spe$slide %in% c("V13B23-339","V13B23-283"), "joint-283-339", spe$slide)
#read in any nnSVG result to get the feature list after filtering
#tmp = read.csv("processed-data/04_feature_selection/per-slide_svgs/V13B23-301_nnSVG-results.csv", row.names=1)
load("processed-data/04_feature_selection/per-sample_spe/V13B23-334_C1.Rdata")
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
df_default = cbind.data.frame(rowData(spe)[names(default),c("gene_id","gene_name")],"dev"=default,
        "rank"=(length(default)+1)-rank(default))

cat("\n\nRunning batch=dummy_slide model...\n")
format(Sys.time(), tz="EST")
batch <- devianceFeatureSelection(rlz, fam="binomial", batch=as.factor(spe$dummy_slide))
df = left_join(df_default, cbind.data.frame(rowData(spe)[rownames(rlz),c("gene_id","gene_name")], "dev"=batch, "rank"=(length(batch)+1)-rank(batch)),
        by=c("gene_id","gene_name"), suffix=c("_default","_batch"))

cat("\n\nSave batch=dummy_slide output...\n")
format(Sys.time(), tz="EST")
write.csv(df, "processed-data/04_feature_selection/test_bindev_batch-dummy-slide.csv",row.names=FALSE)
cat("\nSaved to: processed-data/04_feature_selection/test_bindev_batch-dummy-slide.csv\n")

cat("\n\nRunning batch=sample model...\n")
format(Sys.time(), tz="EST")
batch <- devianceFeatureSelection(rlz, fam="binomial", batch=as.factor(spe$sample_id))
df = left_join(df_default, cbind.data.frame(rowData(spe)[rownames(rlz),c("gene_id","gene_name")], "dev"=batch, "rank"=(length(batch)+1)-rank(batch)),
	by=c("gene_id","gene_name"), suffix=c("_default","_batch")) 

cat("\n\nSave batch=sample output...\n")
format(Sys.time(), tz="EST")
write.csv(df, "processed-data/04_feature_selection/test_bindev_batch-sample.csv",row.names=FALSE)
cat("\nSaved to: processed-data/04_feature_selection/test_bindev_batch-sample.csv\n")

cat("\n\nRunning batch=sex model...\n")
format(Sys.time(), tz="EST")
batch <- devianceFeatureSelection(rlz, fam="binomial", batch=as.factor(spe$sex))
df = left_join(df_default, cbind.data.frame(rowData(spe)[rownames(rlz),c("gene_id","gene_name")], "dev"=batch, "rank"=(length(batch)+1)-rank(batch)),
        by=c("gene_id","gene_name"), suffix=c("_default","_batch")) 

cat("\n\nSave batch=sex output...\n")
format(Sys.time(), tz="EST")
write.csv(df, "processed-data/04_feature_selection/test_bindev_batch-sex.csv",row.names=FALSE)
cat("\nSaved to: processed-data/04_feature_selection/test_bindev_batch-sex.csv\n")

cat("\n\nRunning batch=condition model...\n")
format(Sys.time(), tz="EST")
batch <- devianceFeatureSelection(rlz, fam="binomial", batch=as.factor(spe$condition))
df = left_join(df_default, cbind.data.frame(rowData(spe)[rownames(rlz),c("gene_id","gene_name")], "dev"=batch, "rank"=(length(batch)+1)-rank(batch)),
        by=c("gene_id","gene_name"), suffix=c("_default","_batch")) 

cat("\n\nSave batch=condition output...\n")
format(Sys.time(), tz="EST")
write.csv(df, "processed-data/04_feature_selection/test_bindev_batch-condition.csv",row.names=FALSE)
cat("\nSaved to: processed-data/04_feature_selection/test_bindev_batch-condition.csv\n")


## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
