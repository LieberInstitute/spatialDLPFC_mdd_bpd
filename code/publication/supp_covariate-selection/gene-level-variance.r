setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(pheatmap)
	library(scater)
	library(dplyr)
	library(ggplot2)
})

set.seed(123)

# load precast smoothed
load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
colnames(colData(spe_pseudo))[grep("subsets_mito_percent", colnames(colData(spe_pseudo)))] = "chrM_ratio"
spe_pseudo$pc3 = reducedDim(spe_pseudo, "PCA_1663")[,3]
spe_pseudo$domain = spe_pseudo$smoothed_k9_1663
spe_sm <- spe_pseudo

#need to remove age from colData because of significant digits change
demo = read.csv("processed-data/publication/demographics.csv")
new.cdata = merge(colData(spe_sm)[,setdiff(colnames(colData(spe_sm)), c("age","RIN"))], demo, sort=F)
stopifnot(identical(spe_sm$total, new.cdata$total))
colData(spe_sm) <- new.cdata


# load seurat label
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
colnames(colData(spe_pseudo))[grep("subsets_mito_percent", colnames(colData(spe_pseudo)))] = "chrM_ratio"
spe_pseudo$domain = spe_pseudo$seurat_label
spe_pseudo$pc3 = reducedDim(spe_pseudo, "PCA_1663")[,3]
spe_se <- spe_pseudo
new.cdata = merge(colData(spe_se)[,setdiff(colnames(colData(spe_se)), c("age","RIN"))], demo, sort=F)
stopifnot(identical(spe_se$total, new.cdata$total))
colData(spe_se) <- new.cdata


other.vars = c("sample_id", "nspots", "chrM_ratio", "pc3","age", "BMI", "RIN", "Smoking", "slide", "seq")


var.m = getVarianceExplained(spe_sm, variables=other.vars, exprs_values="logcounts")
summary(var.m)

cor.var.m = cor(var.m, method="pearson")
cor.var.m[cor.var.m==1] = NA

ordered = c("sample_id","slide","pc3","seq","nspots","age","Smoking","BMI","chrM_ratio","RIN")

phm1 = pheatmap(cor.var.m[ordered, ordered], angle_col=90, main="PRECAST (smoothed) gene var corr",
                cluster_cols=F, cluster_rows=F,
                color=colorRampPalette(RColorBrewer::brewer.pal(n=7, "RdBu")[7:1])(100), na_col="grey50",
                breaks=seq(-1, 1, length.out=101))

var.m2 = getVarianceExplained(spe_se, variables=other.vars, exprs_values="logcounts")
summary(var.m2)

cor.var.m2 = cor(var.m2, method="pearson")
cor.var.m2[cor.var.m2==1] = NA

phm2 = pheatmap(cor.var.m2[ordered, ordered], angle_col=90, main="Seurat label gene var corr",
                cluster_cols=F, cluster_rows=F,
                color=colorRampPalette(RColorBrewer::brewer.pal(n=7, "RdBu")[7:1])(100), na_col="grey50",
                breaks=seq(-1, 1, length.out=101))


pdf(file="plots/publication/supp_covariate-selection/gene-level-variance_heatmaps.pdf", width=5, height=5)
plot(phm1[[4]])
plot(phm2[[4]])
dev.off()

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()

