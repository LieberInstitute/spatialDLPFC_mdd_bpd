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

#get sequencing round information
source("code/02_build_spe/getMBvSampleInfo_function.r")
demo = getMBvSampleInfo(REDCapFile="Visium_DATA_2025-01-22_1406.csv",
                        demoFile="DLPFC_cross-disorders_demographics_MBv.csv")
cdata = merge(colData(spe), demo)
#fix "Mbv_034" to "MBv_034"
cdata[grep("b", cdata$MBv_sample),"MBv_sample"] = "MBv_034"
#fix MBv samples with _DO-NO_SEQ
v1 = unique(cdata$MBv_sample)
names(v1) = v1
v2 = sapply(strsplit(v1,"_"), length)
v2[v2>2]
substr(names(v2[v2>2]), start=0, stop=7)
cdata$MBv_sample = substr(cdata$MBv_sample, start=0, stop=7)
#Samples MBv_001-008 and MBv_013-016 were sequenced on an S4 NovaSeq 6000 (Illumina) at the SC-TC. 
#MBv_017-024 were run on a S4 NovaSeq 6000 at the SKCCC. 
#The remaining samples (MBv_009-012 and MBv_025-120) were sequencing at Psomagen over 9 lanes of a 25B Novaseq X. 
#MBv_121-128, and were sequenced on one lane of a 25B Novaseq X at Psomagen in a fourth sequencing batch
cdata$seq = ifelse(cdata$MBv_sample %in% c(paste0("MBv_00",1:8), paste0("MBv_0",13:16)), "SC-TC", "other")
cdata$seq = ifelse(cdata$MBv_sample %in% paste0("MBv_0",17:24), "SKCCC", cdata$seq)
cdata$seq = ifelse(cdata$MBv_sample %in% c("MBv_009", paste0("MBv_0",10:12), paste0("MBv_0",25:99), paste0("MBv_",100:120)), "Psomagen-1", cdata$seq)
cdata$seq = ifelse(cdata$MBv_sample %in% paste0("MBv_",121:128), "Psomagen-2", cdata$seq)
table(cdata[,c("round","seq")])

#reorder cdata to match with spe
rownames(cdata) = cdata$key
cdata = cdata[colnames(spe),]
stopifnot(identical(rownames(cdata), colnames(spe)))
colData(spe) <- cdata


#read in any sample spe and take rownames for subsetting
load("processed-data/04_feature_selection/per-sample_spe/V13B23-334_C1.Rdata")
spe <- spe[rownames(tmp),]
dim(spe)

cat("\nRealize subset spe...\n")
format(Sys.time(), tz="EST")
rlz=writeHDF5Array(counts(spe)) 
getHDF5DumpDir()
format(Sys.time(), tz="EST")
dim(rlz)

#default model
cat("\nRunning default model...\n")
format(Sys.time(), tz="EST")
default <- devianceFeatureSelection(rlz, fam="binomial", batch=NULL)
df_default = cbind.data.frame(rowData(spe)[names(default),c("gene_id","gene_name")],"dev"=default,
        "rank"=(length(default)+1)-rank(default))

#dummy slide (re-put V13B23-283 on 339 slide)
cat("\n\nRunning batch=dummy_slide model...\n")
format(Sys.time(), tz="EST")
batch <- devianceFeatureSelection(rlz, fam="binomial", batch=as.factor(spe$dummy_slide))
df = left_join(df_default, cbind.data.frame(rowData(spe)[rownames(rlz),c("gene_id","gene_name")], "dev"=batch, "rank"=(length(batch)+1)-rank(batch)),
        by=c("gene_id","gene_name"), suffix=c("_default","_batch"))

cat("\n\nSave batch=dummy_slide output...\n")
format(Sys.time(), tz="EST")
write.csv(df, "processed-data/04_feature_selection/test_bindev_batch-dummy-slide.csv",row.names=FALSE)
cat("\nSaved to: processed-data/04_feature_selection/test_bindev_batch-dummy-slide.csv\n")

#sample
cat("\n\nRunning batch=sample model...\n")
format(Sys.time(), tz="EST")
batch <- devianceFeatureSelection(rlz, fam="binomial", batch=as.factor(spe$sample_id))
df = left_join(df_default, cbind.data.frame(rowData(spe)[rownames(rlz),c("gene_id","gene_name")], "dev"=batch, "rank"=(length(batch)+1)-rank(batch)),
	by=c("gene_id","gene_name"), suffix=c("_default","_batch")) 

cat("\n\nSave batch=sample output...\n")
format(Sys.time(), tz="EST")
write.csv(df, "processed-data/04_feature_selection/test_bindev_batch-sample.csv",row.names=FALSE)
cat("\nSaved to: processed-data/04_feature_selection/test_bindev_batch-sample.csv\n")

#seq round
cat("\n\nRunning batch=sequencing round...\n")
format(Sys.time(), tz="EST")
batch <- devianceFeatureSelection(rlz, fam="binomial", batch=as.factor(spe$seq))
df = left_join(df_default, cbind.data.frame(rowData(spe)[rownames(rlz),c("gene_id","gene_name")], "dev"=batch, "rank"=(length(batch)+1)-rank(batch)),
        by=c("gene_id","gene_name"), suffix=c("_default","_batch"))

cat("\n\nSave batch=sequencing round output...\n")
format(Sys.time(), tz="EST")
write.csv(df, "processed-data/04_feature_selection/test_bindev_batch-seq-rnd.csv",row.names=FALSE)
cat("\nSaved to: processed-data/04_feature_selection/test_bindev_batch-seq-rnd.csv\n")

#sex
cat("\n\nRunning batch=sex model...\n")
format(Sys.time(), tz="EST")
batch <- devianceFeatureSelection(rlz, fam="binomial", batch=as.factor(spe$sex))
df = left_join(df_default, cbind.data.frame(rowData(spe)[rownames(rlz),c("gene_id","gene_name")], "dev"=batch, "rank"=(length(batch)+1)-rank(batch)),
        by=c("gene_id","gene_name"), suffix=c("_default","_batch")) 

cat("\n\nSave batch=sex output...\n")
format(Sys.time(), tz="EST")
write.csv(df, "processed-data/04_feature_selection/test_bindev_batch-sex.csv",row.names=FALSE)
cat("\nSaved to: processed-data/04_feature_selection/test_bindev_batch-sex.csv\n")

#condition
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
