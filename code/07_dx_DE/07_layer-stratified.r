setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(DelayedArray)
	#library(spatialLIBD)
	library(edgeR)
	library(BiocParallel)
})

set.seed(123)


load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
#load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
dim(spe_pseudo)
spe_pseudo$pc3 = reducedDim(spe_pseudo)[,"PC3"]

#revised genes, recalculated on pb-filtered samples but not on individual cluster (so that all models from same annotation have same genes tested)
### I checked that the gene filters produced by calculating on the cluster-specific spe object is about the same as this approach (~12k-13k genes)
rowData(spe_pseudo)$high_expr_group_sample_id2 <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)

rowData(spe_pseudo)$high_expr_group_cluster2 <- filterByExpr(spe_pseudo, group = spe_pseudo$smoothed_k9_1663)
#rowData(spe_pseudo)$high_expr_group_cluster2 <- filterByExpr(spe_pseudo, group = spe_pseudo$seurat_label)

spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_sample_id2==T & rowData(spe_pseudo)$high_expr_group_cluster2==T,]
dim(spe_pseudo) 

clust = levels(spe_pseudo$smoothed_k9_1663)
#clust = levels(spe_pseudo$seurat_label)
names(clust) <- clust

fitList <- bplapply(clust, function(x) {
	tmp <- spe_pseudo[,spe_pseudo$smoothed_k9_1663==x]
	#tmp <- spe_pseudo[,spe_pseudo$seurat_label==x]
	#dim(spe_pseudo)

	#make DGE
	dge_pseudo = DGEList(counts(tmp))
	dge_pseudo <- calcNormFactors(dge_pseudo)

	#cat("\ndx model: ~ 0 + group + pc3\n")

	group = interaction(tmp$condition, tmp$sex)
	#table(group)
	dx_mod <- model.matrix(
	  ~ 0 + group + pc3,
	  colData(tmp)
	)
	stopifnot(is.fullrank(dx_mod))

	#cat("\nvoom applied = TRUE\n")
	y = voom(dge_pseudo, dx_mod, plot=F)
	#corfit <- duplicateCorrelation(y, block = colData(spe_pseudo)$sample_id)
	fit <- lmFit(y) #, block = colData(spe_pseudo)$sample_id, correlation = corfit$consensus)
	return(fit)
}, BPPARAM=MulticoreParam(workers=4, RNGseed=1000))

saveRDS(fitList, "processed-data/07_dx_DE/lmFit-voom_layer-stratified_smoothed-k9-1663_condition-sex_rev-gene-input_covars-pc3.rda")
cat("\nlmFit results/ list object saved to: processed-data/07_dx_DE/lmFit-voom_layer-stratified_smoothed-k9-1663_condition-sex_rev-gene-input_covars-pc3.rda\n")

#saveRDS(fitList, "processed-data/07_dx_DE/lmFit-voom_layer-stratified_seurat-pc30_condition-sex_rev-gene-input_covars-pc3.rda")
#cat("\nlmFit results/ list object saved to: processed-data/07_dx_DE/lmFit-voom_layer-stratified_seurat-pc30_condition-sex_rev-gene-input_covars-pc3.rda\n")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()

