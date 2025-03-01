setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(Seurat)
	library(PRECAST)
	library(dplyr)
})
set.seed(123)

load("processed-data/05_clustering/srt-list_spe-HDF5_n1053_counts.Rdata")
#load("processed-data/05_clustering/srt-list_spe-HDF5_n1721_counts.Rdata")
#load("processed-data/05_clustering/srt-list_spe-HDF5_n3198_counts.Rdata")

#use customGeneList to remove 3 RPS|RPL genes present
cglist = setdiff(rownames(srt.sets[[1]]), c("ENSG00000266472","ENSG00000062582","ENSG00000182154"))

#run precast
preobj <- CreatePRECASTObject(seuList = srt.sets, customGenelist=cglist, #rownames(srt.sets[[1]]),
	premin.spots=0, premin.features=0, postmin.spots=0, postmin.features=0)
PRECASTObj <- AddAdjList(preobj, platform = "Visium")
# define model parameters and run model
PRECASTObj <- AddParSetting(PRECASTObj, maxIter = 20, verbose = TRUE, Sigma_equal=FALSE, coreNum=12)
cat("\n\nStart model:",format(Sys.time(), tz="EST"),"\n")
PRECASTObj <- PRECAST(PRECASTObj, K=9)
# pick model (necessary but only changes things if more than 1 K) and integrate
PRECASTObj <- SelectModel(PRECASTObj, criteria="MBIC")
seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")

#save(PRECASTObj, file=here("processed-data","05_clustering","srt_precast-list_k-8-max-it-50_24samp_bindev.3k.Rdata"))
save(seuInt,file="processed-data/05_clustering/srt_precast_120samp_k-9_n1050-repeat.Rdata")

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
