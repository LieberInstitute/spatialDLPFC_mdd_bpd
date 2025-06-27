setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(edgeR)
})
set.seed(123)
#setAutoBlockSize(1e9)

#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
#load("processed-data/06_pseudobulk/spe_n119_pseudo_sample-n1663-k9_norm-filt-combat-seq.Rdata")
load("processed-data/06_pseudobulk/spe_n119_pseudo_sample-n1663-k9_norm-filt.Rdata")
dim(spe_pseudo)

spe_pseudo$pc3 = reducedDim(spe_pseudo)[,"PC3"]

#s.cdata = read.csv("processed-data/06_pseudobulk/continuous_batch_variable_slide-sample-id.csv")
#colData(spe_pseudo) <- merge(colData(spe_pseudo), s.cdata[,c("sample_id","m40_fitted")])

dge_pseudo = DGEList(counts(spe_pseudo))
dge_pseudo <- calcNormFactors(dge_pseudo)

### the way that makes more sense to me for setting up contrasts
dx = spe_pseudo$condition
clus = spe_pseudo$precast_k9_1663
sex = spe_pseudo$sex
group = interaction(dx, clus, sex)
#group = interaction(dx, clus)
cat("\ngroup var produced by interaction():\n")
table(group)

cat("\ndx model: ~ 0 + group + pc3\n")

dx_mod <- model.matrix(
  ~ 0 + group + pc3,
  colData(spe_pseudo)
)
stopifnot(is.fullrank(dx_mod))

#cat("\nvoomLmFit applied = TRUE\n")
#fit <- voomLmFit(dge_pseudo, design=dx_mod, block=spe_pseudo$sample_id)

cat("\nvoom applied = TRUE\n")
y = voom(dge_pseudo, dx_mod, plot=F)
corfit <- duplicateCorrelation(y, block = colData(spe_pseudo)$sample_id)
fit <- lmFit(y, block = colData(spe_pseudo)$sample_id, correlation = corfit$consensus)

#cat("\nvoom applied = FALSE\n")
#corfit <- duplicateCorrelation(logcounts(spe_pseudo), design=dx_mod, block=spe_pseudo$sample_id)
#fit <- lmFit(logcounts(spe_pseudo), design=dx_mod, block=spe_pseudo$sample_id, correlation= corfit$consensus) 

saveRDS(fit, "processed-data/07_dx_DE/lmFit-voom_layer-restricted_precast-k9-1663_condition-sex_covars-pc3.rda")
cat("\nSaved to: processed-data/07_dx_DE/lmFit-voom_layer-restricted_precast-k9-1663_condition-sex_covars-pc3.rda\n")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
