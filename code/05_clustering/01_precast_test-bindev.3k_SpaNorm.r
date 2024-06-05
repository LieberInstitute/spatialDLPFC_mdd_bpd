setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(Seurat)
	library(PRECAST)
	library(here)
})
set.seed(123)

#load srt list and feature list(s)
load(here("processed-data","05_clustering","srt-list_spe_counts.Rdata"))
feature.list = readRDS(here("processed-data","04_preprocessing","bindev-2k-3k_svg_feature-list.rda"))
#load SpaNorm object to access logcounts
load(here("processed-data","04_preprocessing","spe_SpaNorm-test6-bindev.3k.Rdata"))
l1 = unique(spe_small$sample_id)
names(l1) = lapply(l1, function(x) unique(colData(spe_small)[spe_small$sample_id==x,"brain"]))
l1 = lapply(l1, function(x) spe_small[,colData(spe_small)$sample_id==x])
#create preobj with the name elements as the list form of SpaNorm
preobj <- CreatePRECASTObject(seuList = srt.sets[names(l1)], customGenelist=feature.list$bindev.3k, 
	premin.spots=0, premin.features=0, postmin.spots=0, postmin.features=0)
#replace RNA@data with SpaNorm logcounts
for(i in seq_along(l1)) {
	rnames = rownames(preobj@seulist[[i]]$RNA@counts)
	repl.data = logcounts(l1[[i]])[rnames,]
	colnames(repl.data) = paste(l1[[i]]$sample_id, colnames(l1[[i]]), sep="_")
	preobj@seulist[[i]]$RNA@data = repl.data
}

PRECASTObj <- AddAdjList(preobj, platform = "Visium")
# define model parameters and run model
PRECASTObj <- AddParSetting(PRECASTObj, maxIter = 20, verbose = TRUE)
PRECASTObj <- PRECAST(PRECASTObj, K = 7)
# pick model (necessary but only changes things if more than 1 K) and integrate
PRECASTObj <- SelectModel(PRECASTObj, criteria="MBIC")
seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")

save(seuInt,file=here("processed-data","05_clustering","srt_precast_6samp_bindev.3k_SpaNorm.Rdata"))

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
