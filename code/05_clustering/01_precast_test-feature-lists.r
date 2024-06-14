setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(PRECAST)
	library(Seurat)
	library(parallel)
	library(here)
})
set.seed(123)

if("srt-list_spe_counts.Rdata" %in% list.files(here("processed-data","05_clustering"))) {
	cat("\nLoading saved srt.sets object...\n")
	load(here("processed-data","05_clustering","srt-list_spe_counts.Rdata"))
} else {
	load(here("processed-data","04_preprocessing","spe_norm.Rdata"))
	
	l2 = unique(spe$sample_id)
	names(l2) = lapply(l2, function(x) unique(colData(spe)[spe$sample_id==x,"brain"]))
	l2 = lapply(l2, function(x) spe[,colData(spe)$sample_id==x])

	srt.sets = lapply(l2, function(x) {
		rownames(colData(x)) <- paste(x$sample_id, rownames(colData(x)), sep="_")
		colnames(counts(x)) <- rownames(colData(x))
		colData(x)$col <- x$array_col
		colData(x)$row <- x$array_row
		count <- counts(x)
		a1 <- CreateAssayObject(count, assay = "RNA", min.features = 0, min.cells = 0)
		CreateSeuratObject(a1, meta.data = as.data.frame(colData(x)))
	})
	save(srt.sets, file=here("processed-data","05_clustering","srt-list_spe_counts.Rdata"))
}

feature.list = readRDS(here("processed-data","04_preprocessing","bindev-2k-3k_svg_feature-list.rda"))

# for testing purposes
srt.sets = srt.sets[c(1,5,9,13,17,21)]

mclapply(seq_along(feature.list), function(x) {
	cat("\n",names(feature.list)[x],format(Sys.time(),tz="UTC"),"\n")
	preobj <- CreatePRECASTObject(seuList = srt.sets, customGenelist=feature.list[[x]],
		premin.spots=0, premin.features=0, postmin.spots=0, postmin.features=0)
	PRECASTObj <- AddAdjList(preobj, platform = "Visium") 
	PRECASTObj <- AddParSetting(PRECASTObj, maxIter = 20, verbose = TRUE, Sigma_equal=FALSE)
	PRECASTObj <- PRECAST(PRECASTObj, K = 7)
	PRECASTObj <- SelectModel(PRECASTObj, criteria="MBIC")
	seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")
	save(seuInt,file=here("processed-data","05_clustering",paste0("srt_precast_6samp_",names(feature.list)[x],".Rdata")))
}, mc.cores=3)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
