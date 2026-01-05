args = commandArgs(TRUE)
cat("\nCovariate for sensitivity analysis:", args[[1]],"\n")

setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
        library(SpatialExperiment)
        library(DelayedArray)
        library(edgeR)
	library(dplyr)
})

set.seed(123)

load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
#load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
dim(spe_pseudo)
spe_pseudo$pc3 = reducedDim(spe_pseudo)[,"PC3"]

#add extra covars
cdata = read.csv("raw-data/sample_info/Page_BMI_Smoking_CrossDisorder_120825.csv")
sdata = distinct(as.data.frame(colData(spe_pseudo)[,c("brnum","sample_id","condition","sex","age","PMI","RIN")]))
sdata = left_join(sdata, cdata[,c("BrNum","BMI","Smoking")], by=c("brnum"="BrNum"))
new.cdata = merge(colData(spe_pseudo), sdata, sort=F)
stopifnot(identical(spe_pseudo$total, new.cdata$total))
colData(spe_pseudo) <- new.cdata

#revised genes, recalculated on filtered samples
rowData(spe_pseudo)$high_expr_group_sample_id2 <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)

rowData(spe_pseudo)$high_expr_group_cluster2 <- filterByExpr(spe_pseudo, group = spe_pseudo$smoothed_k9_1663)
#rowData(spe_pseudo)$high_expr_group_cluster2 <- filterByExpr(spe_pseudo, group = spe_pseudo$seurat_label)

spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_sample_id2==T & rowData(spe_pseudo)$high_expr_group_cluster2==T,]
dim(spe_pseudo)


#make DGE
dge_pseudo = DGEList(counts(spe_pseudo))
dge_pseudo <- calcNormFactors(dge_pseudo)

#establish model
cat("\nLayer-agnostic model\ndx model: ~ 0 + group + pc3 +", args[[1]], "\n")

group = interaction(spe_pseudo$condition, spe_pseudo$sex)
table(group)
#small work around to use args variable as covariate
colnames(colData(spe_pseudo))[grep(args[[1]], colnames(colData(spe_pseudo)))] <- "test_covar"
dx_mod <- model.matrix(
  ~ 0 + group + pc3 + test_covar,
  colData(spe_pseudo)
)
colnames(dx_mod)[ncol(dx_mod)] = args[[1]]
stopifnot(is.fullrank(dx_mod))

cat("\nvoom applied = TRUE\n")
y = voom(dge_pseudo, dx_mod, plot=F)
corfit <- duplicateCorrelation(y, block = colData(spe_pseudo)$sample_id)
fit <- lmFit(y, block = colData(spe_pseudo)$sample_id, correlation = corfit$consensus)

saveRDS(fit, paste0("processed-data/07_covariate_sensitivity/lmFit-voom_layer-agnostic_smoothed-k9-1663_condition-sex_rev-gene-input_covars-pc3-",
	args[[1]], ".rda"))
cat("\nlmFit results/ object saved to:", paste0("processed-data/07_covariate_sensitivity/lmFit-voom_layer-agnostic_smoothed-k9-1663_condition-sex_rev-gene-input_covars-pc3-",
	args[[1]], ".rda),"\n")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
