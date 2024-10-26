setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
        library(SpatialExperiment)
        library(scry)
        library(dplyr)
        library(here)
})
set.seed(123)

load(here("processed-data","04_preprocessing","spe_norm.Rdata"))

cat("\nRunning default model...\n")
default <- devianceFeatureSelection(spe, fam="binomial", batch=NULL)

df = cbind.data.frame("gene"=rownames(default),"gene_name"=rowData(default)$gene_name,
	"dev"= rowData(default)$binomial_deviance,
	"rank"=(nrow(default)+1)-rank(rowData(default)$binomial_deviance))
        
cat("\nRunning batch model...\n")
batch <- devianceFeatureSelection(spe, fam="binomial", batch=as.factor(spe$slide))

df = left_join(df, cbind.data.frame("gene"=rownames(rowData(batch)), "dev"=rowData(batch)$binomial_deviance,
	"rank"=(nrow(batch)+1)-rank(rowData(batch)$binomial_deviance)),
	by="gene", suffix=c("_default","_slide"))

write.csv(df, here("processed-data","04_preprocessing","bindev_default-slide.csv"),row.names=FALSE)

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
