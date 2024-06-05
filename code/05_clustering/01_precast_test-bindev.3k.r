setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(Seurat)
	library(PRECAST)
	library(here)
})
set.seed(123)

load(here("processed-data","05_clustering","srt-list_spe_counts.Rdata"))
feature.list = readRDS(here("processed-data","04_preprocessing","bindev-2k-3k_svg_feature-list.rda"))

preobj <- CreatePRECASTObject(seuList = srt.sets, customGenelist=feature.list$bindev.3k,
	premin.spots=0, premin.features=0, postmin.spots=0, postmin.features=0)
PRECASTObj <- AddAdjList(preobj, platform = "Visium")
# define model parameters and run model
PRECASTObj <- AddParSetting(PRECASTObj, maxIter = 20, verbose = TRUE)
PRECASTObj <- PRECAST(PRECASTObj, K = 7)
# pick model (necessary but only changes things if more than 1 K) and integrate
PRECASTObj <- SelectModel(PRECASTObj, criteria="MBIC")
seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")

save(seuInt,file=here("processed-data","05_clustering","srt_precast_24samp_bindev.3k.Rdata"))

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
