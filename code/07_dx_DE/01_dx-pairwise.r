setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(DelayedArray)
	library(spatialLIBD)
})

set.seed(123)

load("processed-data/06_pseudobulk/spe_n119_pseudo_sample-comb-clus_norm-filt.Rdata")
dim(spe_pseudo) # 12397   950
spe_pseudo$slide_name = gsub("-","\\.", spe_pseudo$slide)

dx_mod <- registration_model(spe_pseudo,
    covars = c("detected", "ncells", "sex", "combined_cluster", "age"),
    var_registration = "condition"
)

dx_block_cor <- registration_block_cor(spe_pseudo, registration_model = dx_mod,
    var_sample_id = "sample_id"
)

dx_res <- registration_stats_pairwise(spe_pseudo, block_cor = dx_block_cor,
  #covars = c("age", "sex", "combined_cluster"),
  registration_model= dx_mod,
  var_registration = "condition",
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)

write.csv(dx_res, "processed-data/07_pseudobulk/results_dx-pairwise_covars-age-cluster-detected-ncells-sex.csv", row.names=T)


cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
