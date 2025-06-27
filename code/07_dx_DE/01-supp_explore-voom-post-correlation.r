#run in interactive mode
library(SpatialExperiment)
library(edgeR)
library(dplyr)
library(ggplot2)

set.seed(123)
load("processed-data/06_pseudobulk/spe_n119_pseudo_sample-comb-clus_norm-filt.Rdata")


layer_mod <- spatialLIBD::registration_model(spe_pseudo,
       covars = c("sex","age"),
       var_registration = "combined_cluster"
)
layer_mod

#duplivateCorrelation then voom then duplicate correlation then voom
dge_pseudo = DGEList(counts(spe_pseudo))
dge_pseudo <- calcNormFactors(dge_pseudo)

dx_mod <- model.matrix(
  ~ 0 + condition + combined_cluster + sex,
  colData(spe_pseudo)
)
par(mfrow=c(2,1))
y = voom(dge_pseudo, dx_mod, plot=T)
str(y)
y$weights[1:3,1:4]
#read in correlation from already run model
fit1 <- readRDS("processed-data/07_dx_DE/lmFit-voom_layer-adjusted_condition_covars-cluster-sex.rda")
fit1$correlation #0.3800061

y1 = voom(dge_pseudo, dx_mod, block=spe_pseudo$sample_id, correlation=fit1$correlation, plot=T)
#can see substantially reduced scatter/noise
y1$weights[1:3,1:4]
#can see that weights are changed
