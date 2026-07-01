setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(pheatmap)
	library(scater)
	library(dplyr)
	library(ggplot2)
})

set.seed(123)

#need to remove age from colData because of significant digits change
demo = read.csv("processed-data/publication/supp_tables/demographics.csv")

# collinearity of donor-level variables
cat("\nCorrelation of donor-level variables...\n")
cor1 = cor(demo[,c("age","RIN","BMI")])
cor1[cor1==1] = NA

print(signif(cor1, digits=3))

# load precast smoothed
cat("\nDomain-SP covariate correlations...\n")
load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
colnames(colData(spe_pseudo))[grep("subsets_mito_percent", colnames(colData(spe_pseudo)))] = "chrM_ratio"
spe_pseudo$pc3 = reducedDim(spe_pseudo, "PCA_1663")[,3]
spe_sm <- spe_pseudo

# collinearity of pseudobulk-level variables
cor2 = cor(as.data.frame(colData(spe_sm)[,c("sum","detected","nspots","chrM_ratio")]))
cor2[cor2==1] = NA

cat("Pseudobulk level variables...\n")
print(signif(cor2, digits=3))

new.cdata = merge(colData(spe_sm)[,c("sample_id","brnum","condition","sex","nspots","chrM_ratio","pc3")], 
	demo[,c("sample_id","brnum","sex","age","RIN","BMI","Smoking")], sort=F)
stopifnot(identical(spe_sm$pc3, new.cdata$pc3))

cor3 = cor(as.data.frame(new.cdata[,c("age","RIN","BMI","pc3","nspots","chrM_ratio")]))
cor3[cor3==1] = NA

cat("\nFinal covariate selection...\n")
print(signif(cor3, digits=3))


# load seurat label
cat("\nDomain-CT covariate correlations...\n")
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
colnames(colData(spe_pseudo))[grep("subsets_mito_percent", colnames(colData(spe_pseudo)))] = "chrM_ratio"
spe_pseudo$pc3 = reducedDim(spe_pseudo, "PCA_1663")[,3]
spe_se <- spe_pseudo

cor2 = cor(as.data.frame(colData(spe_se)[,c("sum","detected","nspots","chrM_ratio")]))
cor2[cor2==1] = NA

cat("Pseudobulk	level variables...\n")
print(signif(cor2, digits=3))

new.cdata = merge(colData(spe_se)[,c("sample_id","brnum","condition","sex","nspots","chrM_ratio","pc3")], 
        demo[,c("sample_id","brnum","sex","age","RIN","BMI","Smoking")], sort=F)
stopifnot(identical(spe_se$pc3, new.cdata$pc3))

cor3 = cor(as.data.frame(new.cdata[,c("age","RIN","BMI","pc3","nspots","chrM_ratio")]))
cor3[cor3==1] = NA

cat("\nFinal covariate selection...\n")
print(signif(cor3, digits=3))


cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
