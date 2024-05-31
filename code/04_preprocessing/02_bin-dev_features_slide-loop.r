setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(scry)
	library(dplyr)
	library(parallel)
	library(here)
})
set.seed(123)

load(here("processed-data","04_preprocessing","spe_norm.Rdata"))

l1 = unique(spe$slide)
names(l1) = l1
l1 = lapply(l1, function(x) spe[-grep("^MT-",rowData(spe)$gene_name), spe$slide==x])

mclapply(l1, function(x) {
	cat("\n",unique(x$slide),": Running default model...\n")
        default <- devianceFeatureSelection(x, fam="binomial", batch=NULL)

        df = cbind.data.frame("gene"=rownames(default),"gene_name"=rowData(default)$gene_name,
                "dev"= rowData(default)$binomial_deviance,
                "rank"=(nrow(default)+1)-rank(rowData(default)$binomial_deviance))
	
	cat("\n",unique(x$slide),": Running batch model...\n")
        batch <- devianceFeatureSelection(x, fam="binomial", batch=as.factor(x$brain))

        df = left_join(df, cbind.data.frame("gene"=rownames(rowData(batch)), "dev"=rowData(batch)$binomial_deviance,
                "rank"=(nrow(batch)+1)-rank(rowData(batch)$binomial_deviance)),
        by="gene", suffix=c("_default","_brain"))

        write.csv(df, here("processed-data","04_preprocessing",paste0("bindev_",unique(x$slide),"_default-brain.csv")),row.names=FALSE)
}, mc.cores=6)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
