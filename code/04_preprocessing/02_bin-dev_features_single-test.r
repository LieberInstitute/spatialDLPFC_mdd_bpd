setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(scry)
	library(scran)
	library(dplyr)
	library(here)
})
set.seed(123)

load(here("processed-data","04_preprocessing","spe_norm.Rdata"))
table(spe$plot_outliers)
l1 = unique(spe$slide)

#pick 4 samples (1 slide), 8 samples (2 slides), or 16 samples (4slides)
spe_small = spe[,spe$slide %in% l1[c(1,3)]]
spe_small = spe_small[rowSums(logcounts(spe_small))>3,]
dim(spe_small)
table(spe_small$brain)

#first run without correction
cat("\nExecuting default model... \n")
dev <- nullResiduals(counts(spe_small), type="deviance", fam="binomial")
#dev = assays(spe_small)$binomial_deviance_residuals

df = cbind.data.frame("gene"=rownames(dev),"gene_name"=rowData(spe_small)$gene_name,
	"dev"= matrixStats::rowVars(as.matrix(dev)),
	"rank"=(nrow(dev)+1)-rank(matrixStats::rowVars(as.matrix(dev))))
#then run with correction
cat("\nExecuting batch model... \n")
batch_test <- devianceFeatureSelection(spe_small, fam="binomial", batch=as.factor(spe_small$brain))

df = left_join(df, cbind.data.frame("gene"=rownames(rowData(batch_test)),
	"dev"=rowData(batch_test)$binomial_deviance,
	"rank"=(nrow(batch_test)+1)-rank(rowData(batch_test)$binomial_deviance)),
	by="gene", suffix=c("_default","_brain"))

write.csv(df, here("processed-data","04_preprocessing","bindev_test-single_n8_default-brain.csv"), row.names=FALSE)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
