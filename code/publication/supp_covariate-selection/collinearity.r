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
spe_sm <- spe_pseudo

#need to remove age from colData because of significant digits change
demo = read.csv("processed-data/publication/demographics.csv")
new.cdata = merge(colData(spe_sm)[,setdiff(colnames(colData(spe_sm)), c("age","RIN"))], demo, sort=F)
stopifnot(identical(spe_sm$total, new.cdata$total))
colData(spe_sm) <- new.cdata


# load seurat label
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
colnames(colData(spe_pseudo))[grep("subsets_mito_percent", colnames(colData(spe_pseudo)))] = "chrM_ratio"
spe_pseudo$pc3 = reducedDim(spe_pseudo, "PCA_1663")[,3]
spe_se <- spe_pseudo
new.cdata = merge(colData(spe_se)[,setdiff(colnames(colData(spe_se)), c("age","RIN"))], demo, sort=F)
stopifnot(identical(spe_se$total, new.cdata$total))
colData(spe_se) <- new.cdata
other.vars = c("nspots", "chrM_ratio", "pc3","age", "BMI", "RIN", "Smoking", "slide", "seq")

# collinearity of donor-level variables
cor1 = cor(demo[,c("age","RIN","BMI")])
cor1_num = signif(cor1, digits=3)
cor1[cor1==1] = NA

phm1 = pheatmap(cor1, display_numbers= cor1_num, cluster_rows=F, cluster_cols=F, main="Donor-level variables",
         color=colorRampPalette(RColorBrewer::brewer.pal(n=7, "RdBu")[7:1])(100),
         breaks=seq(-1, 1, length.out=101), na_col = "grey50", angle_col=0)

# collinearity of pseudobulk-level variables
cor2 = cor(as.data.frame(colData(spe_sm)[,c("sum","detected","nspots","chrM_ratio","pc3")]))
cor2[cor2==1] = NA
phm2 = pheatmap(cor2, cluster_rows=F, cluster_cols=F, main="PRECAST (smoothed) pb variables",
         color=colorRampPalette(RColorBrewer::brewer.pal(n=7, "RdBu")[7:1])(100),
         breaks=seq(-1, 1, length.out=101), na_col = "grey50", angle_col=0)

cor3 = cor(as.data.frame(colData(spe_se)[,c("sum","detected","nspots","chrM_ratio","pc3")]))
cor3[cor3==1] = NA
phm3 = pheatmap(cor3, cluster_rows=F, cluster_cols=F, main="Seurat label pb variables",
         color=colorRampPalette(RColorBrewer::brewer.pal(n=7, "RdBu")[7:1])(100),
         breaks=seq(-1, 1, length.out=101), na_col = "grey50", angle_col=0)


pdf(file="plots/publication/supp_covariate-selection/continuous-variable_collinearity.pdf", width=4, height=4)
plot(phm1[[4]])
plot(phm2[[4]])
plot(phm3[[4]])
dev.off()

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
