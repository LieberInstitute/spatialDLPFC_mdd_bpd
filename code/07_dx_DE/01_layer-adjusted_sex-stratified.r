setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(DelayedArray)
	#library(spatialLIBD)
	library(edgeR)
})

set.seed(123)

load("processed-data/06_pseudobulk/spe_n119_pseudo_sample-n1663-k9_norm-filt.Rdata")
#load("processed-data/06_pseudobulk/spe_n119_pseudo_sample-n1663-k9_norm-filt-combat-seq.Rdata")
dim(spe_pseudo)

s.cdata = read.csv("processed-data/06_pseudobulk/continuous_batch_variable_slide-sample-id.csv")
colData(spe_pseudo) <- merge(colData(spe_pseudo), s.cdata[,c("sample_id","m40_fitted","k3")])
spe_pseudo$k3 = as.factor(spe_pseudo$k3)

#stratify by sex
spe_f = spe_pseudo[,spe_pseudo$sex=="F"]
spe_m = spe_pseudo[,spe_pseudo$sex=="M"]

spe_list = list(spe_f, spe_m)


#analyze both in for loop (so that when testing different covars only have to change single dx_mod) 
for(i in spe_list) {

spe_tmp = i
test_sex = unique(spe_tmp$sex)

cat("\n\nANALYZING SEX =",test_sex,"\n")
dim(spe_tmp)

dge_pseudo = DGEList(counts(spe_tmp))
dge_pseudo <- calcNormFactors(dge_pseudo)

cat("\ndx model: ~ 0 + condition + precast_k9_1663 + m40_fitted\n")

dx_mod <- model.matrix(
  ~ 0 + condition + precast_k9_1663 + m40_fitted,
  colData(spe_tmp)
)
stopifnot(is.fullrank(dx_mod))

#cat("\nvoomLmFit applied = TRUE\n")
#fit <- voomLmFit(dge_pseudo, design=dx_mod, block=spe_pseudo$sample_id)

cat("\nvoom applied = TRUE\n")
y = voom(dge_pseudo, dx_mod, plot=F)
corfit <- duplicateCorrelation(y, block = colData(spe_tmp)$sample_id)
fit <- lmFit(y, block = colData(spe_tmp)$sample_id, correlation = corfit$consensus)

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

saveRDS(fit, paste0("processed-data/07_dx_DE/lmFit-voom_layer-adjusted_precast-k9-1663_condition_", test_sex, "-only_covars-cluster-m40.rda"))
cat("\nlmFit results/ object saved to:", 
	paste0("processed-data/07_dx_DE/lmFit-voom_layer-adjusted_precast-k9-1663_condition_", 
		test_sex, "-only_covars-cluster-m40.rda"),
	"\n")

}

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
