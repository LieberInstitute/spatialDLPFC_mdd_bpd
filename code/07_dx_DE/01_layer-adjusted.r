setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(DelayedArray)
	#library(spatialLIBD)
	library(edgeR)
})

set.seed(123)

#load("processed-data/06_pseudobulk/spe_n119_pseudo_sample-comb-clus_norm-filt.Rdata")
load("processed-data/06_pseudobulk/spe_n119_pseudo_sample-n1663-k9_norm-filt.Rdata")
#load("processed-data/06_pseudobulk/spe_n119_pseudo_sample-n1663-k9_norm-filt-combat-seq.Rdata")
dim(spe_pseudo) # 12397   950
spe_pseudo$slide_name = gsub("-","\\.", spe_pseudo$slide)

##vascular integrity score
#bbb.stroke=c("HAMP","CCL2","ZFP36","C11orf96","GADD45B")
#angio = c("ANGPTL4","VEGFA","CHI3L1","SERPINA3")
#vasc.comb = c(angio, bbb.stroke, "MT1X")
#
#vasc.comb.id = rownames(spe_pseudo)[rowData(spe_pseudo)$gene_name %in% vasc.comb]
#colData(spe_pseudo)$vasc.score = scale(colSums(logcounts(spe_pseudo)[vasc.comb.id,]))
#
#bbb.stroke2 = c("HAMP","CCL2","ZFP36","GADD45B")
#bbb.id = rownames(spe_pseudo)[rowData(spe_pseudo)$gene_name %in% bbb.stroke2]
#colData(spe_pseudo)$bbb.score = scale(colSums(logcounts(spe_pseudo)[bbb.id,]))

dge_pseudo = DGEList(counts(spe_pseudo))
dge_pseudo <- calcNormFactors(dge_pseudo)

cat("\ndx model: ~ 0 + condition + precast_k9_1663 + sex + slide_name\n")

dx_mod <- model.matrix(
  ~ 0 + condition + precast_k9_1663 + sex + slide_name,
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

saveRDS(fit, "processed-data/07_dx_DE/lmFit-voom_layer-adjusted_precast-k9-1663_condition_covars-cluster-sex-slide.rda")
cat("\nlmFit results/ object saved to: processed-data/07_dx_DE/lmFit-voom_layer-adjusted_precast-k9-1663_condition_covars-cluster-sex-slide.rda\n")

#cont_mtx = matrix(0, ncol=3, nrow=ncol(coef(fit)), dimnames=list(colnames(coef(fit)), c("NTC.MDD","NTC.BPD","MDD.BPD")))
#cont_mtx["conditionNTC",c("NTC.MDD","NTC.BPD")] <- c(-1, -1)
#cont_mtx["conditionMDD",c("NTC.MDD","MDD.BPD")] <- c(1, -1)
#cont_mtx["conditionBPD",c("NTC.BPD","MDD.BPD")] <- c(1, 1)
#head(cont_mtx)

#tmp1 = contrasts.fit(fit, cont_mtx)
#tmp2 = eBayes(tmp1)

#saveRDS(tmp2, "processed-data/07_dx_DE/eBayes_fit_layer-adjusted_condition_covars-age-cluster-detected-sex-slide.rda")



cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
