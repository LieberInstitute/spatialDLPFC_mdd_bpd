setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(DelayedArray)
	#library(spatialLIBD)
	library(edgeR)
})

set.seed(123)


#load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
dim(spe_pseudo)
spe_pseudo$pc3 = reducedDim(spe_pseudo)[,"PC3"]

#revised genes, recalculated on filtered samples
rowData(spe_pseudo)$high_expr_group_sample_id2 <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
#rowData(spe_pseudo)$high_expr_group_cluster2 <- filterByExpr(spe_pseudo, group = spe_pseudo$smoothed_k9_1663)
rowData(spe_pseudo)$high_expr_group_cluster2 <- filterByExpr(spe_pseudo, group = spe_pseudo$seurat_label)
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_sample_id2==T & rowData(spe_pseudo)$high_expr_group_cluster2==T,]
dim(spe_pseudo) 

##remove lowly expressed genes that were included for comparison with snRNAseq data
#spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_sample_id==T & rowData(spe_pseudo)$high_expr_group_cluster==T,]
#dim(spe_pseudo)

#make DGE
dge_pseudo = DGEList(counts(spe_pseudo))
dge_pseudo <- calcNormFactors(dge_pseudo)

#cat("\ndx model: ~ 0 + group + smoothed_k9_1663 + pc3 + age\n")
cat("\ndx model: ~ 0 + group + seurat_label + pc3 + age\n")

group = interaction(spe_pseudo$condition, spe_pseudo$sex)
table(group)
dx_mod <- model.matrix(
  #~ 0 + group + smoothed_k9_1663 + pc3 +age,
  ~ 0 + group + seurat_label + pc3 + age,
  colData(spe_pseudo)
)
stopifnot(is.fullrank(dx_mod))

#cat("\nvoomLmFit applied = TRUE\n")
#fit <- voomLmFit(dge_pseudo, design=dx_mod, block=spe_pseudo$sample_id)

cat("\nvoom applied = TRUE\n")
y = voom(dge_pseudo, dx_mod, plot=F)
corfit <- duplicateCorrelation(y, block = colData(spe_pseudo)$sample_id)
fit <- lmFit(y, block = colData(spe_pseudo)$sample_id, correlation = corfit$consensus)

#cat("\nvoom re-applied after duplicateCorrelation = TRUE\n")
## https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf PAGE 127-128
#### "The intra cell line correlation will change the voom weights slightly, so we run voom a second time"
### I looked at the effect on the voom mean-variance plot and re-applying voom after calculating the correlation significantly reduces the noise
#y1 = voom(dge_pseudo, design = dx_mod, block = spe_pseudo$sample_id, correlation = corfit$consensus)
##it wants me to recalculate the correlation but that step takes forever so i'm skipping (also b/c in their example it doesn't change)
#fit <- lmFit(y1, block = colData(spe_pseudo)$sample_id, correlation = corfit$consensus)

#cat("\nvoom applied = FALSE\n")
#corfit <- duplicateCorrelation(logcounts(spe_pseudo), design=dx_mod, block=spe_pseudo$sample_id)
#fit <- lmFit(logcounts(spe_pseudo), design=dx_mod, block=spe_pseudo$sample_id, correlation=corfit$consensus)

#saveRDS(fit, "processed-data/07_dx_DE/lmFit-voom_layer-adjusted_smoothed-k9-1663_condition-sex_rev-gene-input_covars-pc3-age.rda")
#cat("\nlmFit results/ object saved to: processed-data/07_dx_DE/lmFit-voom_layer-adjusted_smoothed-k9-1663_condition-sex_rev-gene-input_covars-pc3-age.rda\n")

saveRDS(fit, "processed-data/07_dx_DE/lmFit-voom_layer-adjusted_seurat-pc30_condition-sex_rev-gene-input_covars-pc3-age.rda")
cat("\nlmFit results/ object saved to: processed-data/07_dx_DE/lmFit-voom_layer-adjusted_seurat-pc30_condition-sex_rev-gene-input_covars-pc3-age.rda\n")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
