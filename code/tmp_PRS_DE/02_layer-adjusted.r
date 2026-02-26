args = commandArgs(TRUE)
cat("\nGWAS set for PRS predictor:", args[[1]],"\n")

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
cdata = read.csv("raw-data/PRS/PRS_chosen-p-cutoffs.csv")
cdata = cdata[,c(1:2, 5:6, grep(args[[1]], colnames(cdata)))]
colnames(cdata)[5] = paste0("prs", args[[1]])
sdata = distinct(as.data.frame(colData(spe_pseudo)[,c("brnum","sample_id","condition","sex","age","PMI","RIN")]))
sdata = left_join(sdata, cdata, by=c("brnum","age","PMI","RIN"))
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
cat("\nLayer-adjusted model\ndx model: ~ 0 + ", paste0("prs", args[[1]]), "*sex + smoothed_k9_1663 + pc3 + age + nspots\n")
#cat("\nLayer-adjusted model\ndx model: ~ 0 + ", paste0("prs", args[[1]]), ":sex + seurat_label + pc3\n")

spe_pseudo$age = scale(spe_pseudo$age)
spe_pseudo$nspots = scale(spe_pseudo$nspots)
#small work around to use args variable as covariate
colnames(colData(spe_pseudo))[grep(args[[1]], colnames(colData(spe_pseudo)))] <- "PRS"
spe_pseudo$PRS = scale(spe_pseudo$PRS)
dx_mod <- model.matrix(
  ~ 0 + PRS*sex + smoothed_k9_1663 + pc3 + age + nspots,
  #~ 0 + PRS:sex + seurat_label + pc3,
  colData(spe_pseudo)
)
colnames(dx_mod) <- gsub("PRS",paste0("prs",args[[1]]), colnames(dx_mod))
stopifnot(is.fullrank(dx_mod))

cat("\nvoom applied = TRUE\n")
y = voom(dge_pseudo, dx_mod, plot=F)
corfit <- duplicateCorrelation(y, block = colData(spe_pseudo)$sample_id)
fit <- lmFit(y, block = colData(spe_pseudo)$sample_id, correlation = corfit$consensus)

saveRDS(fit, paste0("processed-data/09_PRS_DE/lmFit-voom_layer-adjusted_smoothed-k9-1663_", paste0("prs", args[[1]]), 
	"-sex_rev-gene-input_covars-pc3-age-nspots.rda"))
cat("\nlmFit results/ object saved to:", paste0("processed-data/09_PRS_DE/lmFit-voom_layer-adjusted_smoothed-k9-1663_", paste0("prs", args[[1]]),
	"-sex_rev-gene-input_covars-pc3-age-nspots.rda"),"\n")

#saveRDS(fit, paste0("processed-data/09_PRS_DE/lmFit-voom_layer-adjusted_seurat-pc30_", paste0("prs", args[[1]]),
#	"-sex_rev-gene-input_covars-pc3.rda"))
#cat("\nlmFit results/ object saved to:", paste0("processed-data/09_PRS_DE/lmFit-voom_layer-adjusted_seurat-pc30_", paste0("prs", args[[1]]),
#	"-sex_rev-gene-input_covars-pc3.rda),"\n")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
