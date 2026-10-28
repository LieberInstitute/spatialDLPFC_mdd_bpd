setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SingleCellExperiment)
	library(edgeR)
})
set.seed(123)

load("processed-data/08_compare_DEGs/SZBDMulti-seq/sce_SZBD_pseudo_indivID-low-res_norm-filt.Rdata")
dim(sce_pseudo)
sce_pseudo$pc10 = reducedDim(sce_pseudo, "PCA_HVG")[,"PC10"]

#revised genes, recalculated on filtered samples
rowData(sce_pseudo)$high_expr_group_sample2 <- filterByExpr(sce_pseudo, group = sce_pseudo$individualID)
rowData(sce_pseudo)$high_expr_group_cluster2 <- filterByExpr(sce_pseudo, group = sce_pseudo$seurat_low.res)
sce_pseudo <- sce_pseudo[rowData(sce_pseudo)$high_expr_group_sample2==T & rowData(sce_pseudo)$high_expr_group_cluster2==T,]
dim(sce_pseudo)

#make DGE
dge_pseudo = DGEList(counts(sce_pseudo))
dge_pseudo <- calcNormFactors(dge_pseudo)

cat("\ndx model: ~ 0 + group + seurat_low.res + pc10 + age + ncells\n")
group = interaction(sce_pseudo$condition, sce_pseudo$sex)
table(group)

sce_pseudo$ncells = scale(sce_pseudo$ncells)
sce_pseudo$age = scale(sce_pseudo$Age_death)

dx_mod <- model.matrix(
  ~ 0 + group + seurat_low.res + pc10 + age + ncells,
  colData(sce_pseudo)
)
stopifnot(is.fullrank(dx_mod))

cat("\nvoom applied = TRUE\n")
y = voom(dge_pseudo, dx_mod, plot=F)
corfit <- duplicateCorrelation(y, block = colData(sce_pseudo)$individualID)
fit <- lmFit(y, block = colData(sce_pseudo)$individualID, correlation = corfit$consensus)

saveRDS(fit, "processed-data/08_compare_DEGs/SZBDMulti-seq/lmFit-voom_layer-adjusted_low-res_condition-sex_covars-pc10-age-ncells.rda")
cat("\nlmFit results/ object saved to: processed-data/08_compare_DEGs/SZBDMulti-seq/lmFit-voom_layer-adjusted_low-res_condition-sex_covars-pc10-age-ncells.rda\n")


# now make layer-restricted model
dx = sce_pseudo$condition
clus = sce_pseudo$seurat_low.res
sex = sce_pseudo$sex
group = interaction(dx, clus, sex)

cat("\ngroup var produced by interaction():\n")
table(group)

cat("\ndx model: ~ 0 + group + pc10 + age + ncells\n")
dx_mod <- model.matrix(
  ~ 0 + group + pc10 + age + ncells,
  colData(sce_pseudo)
)
stopifnot(is.fullrank(dx_mod))

cat("\nvoom applied = TRUE\n")
y = voom(dge_pseudo, dx_mod, plot=F)
corfit <- duplicateCorrelation(y, block = colData(sce_pseudo)$individualID)
fit <- lmFit(y, block = colData(sce_pseudo)$individualID, correlation = corfit$consensus)

saveRDS(fit, "processed-data/08_compare_DEGs/SZBDMulti-seq/lmFit-voom_layer-restricted_low-res_condition-sex_covars-pc10-age-ncells.rda")
cat("\nlmFit results/ object saved to: processed-data/08_compare_DEGs/SZBDMulti-seq/lmFit-voom_layer-restricted_low-res_condition-sex_covars-pc10-age-ncells.rda\n")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()

