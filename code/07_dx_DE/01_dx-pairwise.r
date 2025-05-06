setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(DelayedArray)
	#library(spatialLIBD)
	library(edgeR)
})

set.seed(123)

load("processed-data/06_pseudobulk/spe_n119_pseudo_sample-comb-clus_norm-filt.Rdata")
dim(spe_pseudo) # 12397   950
spe_pseudo$slide_name = gsub("-","\\.", spe_pseudo$slide)

#dx = spe_pseudo$condition
#sex = spe_pseudo$sex
#group = interaction(dx, sex)

dx_mod <- model.matrix(
  ~ 0 + condition + age + combined_cluster + detected + sex + slide_name,
  colData(spe_pseudo)
)

corfit <- duplicateCorrelation(
  logcounts(spe_pseudo),
  design = dx_mod,
  block = colData(spe_pseudo)$sample_id
)

fit <- lmFit(
  logcounts(spe_pseudo),
  design = dx_mod,
  block = colData(spe_pseudo)$sample_id,
  correlation = corfit$consensus
)

cont_mtx = matrix(0, ncol=3, nrow=ncol(coef(fit)), dimnames=list(colnames(coef(fit)), c("NTC.MDD","NTC.BPD","MDD.BPD")))
cont_mtx["conditionNTC",c("NTC.MDD","NTC.BPD")] <- c(-1, -1)
cont_mtx["conditionMDD",c("NTC.MDD","MDD.BPD")] <- c(1, -1)
cont_mtx["conditionBPD",c("NTC.BPD","MDD.BPD")] <- c(1, 1)
head(cont_mtx)

tmp1 = contrasts.fit(fit, cont_mtx)
tmp2 = eBayes(tmp1)

saveRDS(tmp2, "processed-data/07_dx_DE/eBayes_fit_layer-adjusted_condition_covars-age-cluster-detected-sex-slide.rda")

#dx_mod <- registration_model(spe_pseudo,
#    covars = c("detected", "ncells", "sex", "combined_cluster", "age"),
#    var_registration = "condition"
#)

#dx_block_cor <- registration_block_cor(spe_pseudo, registration_model = dx_mod,
#    var_sample_id = "sample_id"
#)

#dx_res <- registration_stats_pairwise(spe_pseudo, block_cor = dx_block_cor,
#  #covars = c("age", "sex", "combined_cluster"),
#  registration_model= dx_mod,
#  var_registration = "condition",
#  gene_ensembl = "gene_id",
#  gene_name = "gene_name"
#)

#write.csv(dx_res, "processed-data/07_pseudobulk/results_dx-pairwise_covars-age-cluster-detected-ncells-sex.csv", row.names=T)


cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
