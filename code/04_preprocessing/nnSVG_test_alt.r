setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')

library(SpatialExperiment)
library(scran)
library(nnSVG)
library(parallel)
library(here)
set.seed(123)

load(here("processed-data","03_QC","spe_demo-filt.Rdata"))

spe$exclude_outliers = spe$sum_local.outlier | spe$genes_local.outlier | spe$chrM.ratio_local.outlier | spe$exclude_east_edge
spe$extra_3MAD_outliers = spe$sum_3MAD.outlier_sample | spe$genes_3MAD.outlier_sample
spe$plot_outliers = factor(paste(spe$exclude_outliers, spe$extra_3MAD_outliers),
	levels=c("FALSE FALSE","FALSE TRUE","TRUE FALSE","TRUE TRUE"),
	labels=c("keep","3MAD flag","exclude","exclude"))

spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)

l1 = unique(spe$slide)
names(l1) = l1
l1 = lapply(l1, function(x) spe[-grep("^MT-",rowData(spe)$gene_name), spe$slide==x & spe$plot_outliers!="exclude"])

for (i in seq_along(l1)) {
	spe_small <- filter_genes(l1[[i]], filter_genes_ncounts = 3, filter_genes_pcspots = .5, filter_mito=T)
	spe_nnSVG <- nnSVG(spe_small, n_threads=12)
	svg = rowData(spe_nnSVG)[rowData(spe_nnSVG)$padj <= 0.05,]
	saveRDS(svg, paste0("/users/jthompso/nnSVG_test_",names(l1)[i],".rda"))
}
