setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
        library(SpatialExperiment)
        library(DelayedArray)
        library(spatialLIBD)
        library(dplyr)
        library(ggplot2)
        library(pheatmap)
})

set.seed(123)

load("processed-data/06_pseudobulk/spe_n119_pseudo_sample-comb-clus_norm-filt.Rdata")
dim(spe_pseudo) # 12397   942

layer_mod <- registration_model(spe_pseudo,
       covars = c("sex","slide","detected","ncells", "age"),
       var_registration = "combined_cluster"
)

layer_block_cor <- registration_block_cor(spe_pseudo, registration_model = layer_mod,
    var_sample_id = "sample_id"
)

layer_res <- registration_stats_enrichment(spe_pseudo, block_cor = layer_block_cor,
  covars = c("sex","slide","detected","ncells", "age"),
  var_registration = "combined_cluster",
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)

write.csv(layer_res, "processed-data/06_pseudobulk/results_layer-enrichment_covars-age-detected-ncells-sex-slide.csv", row.names=T)
cat("\nLayer enrichment test results saved to: processed-data/06_pseudobulk/results_layer-enrichment_covars-age-detected-ncells-sex-slide.csv\n")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
