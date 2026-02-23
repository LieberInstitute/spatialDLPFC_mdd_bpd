setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(dplyr)
})

set.seed(123)

load("processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo-dotplot_azimuth-super-broad.Rdata")
source("code/06_pseudobulk/SZBDMulti-seq/get-mean-ratio_mod.r")
mratio.sn <- get_mean_ratio_mod(sce_summ, "azimuth_super.broad", min_prop.detected= .1,
                           gene_ensembl = "gene_id", gene_name = "gene_name")

write.csv(mratio.sn, "processed-data/06_pseudobulk/SZBDMulti-seq/SZBD-control_azimuth-super-broad_mean-ratio.csv", row.names=F)

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
