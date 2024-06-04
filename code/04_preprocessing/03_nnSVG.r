setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')

suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(scran)
	library(nnSVG)
	library(parallel)
	library(here)
})
set.seed(123)

load(here("processed-data","04_preprocessing","spe_norm.Rdata"))
spe = filter_genes(spe, filter_genes_ncounts = 3, filter_genes_pcspots = .5, filter_mito=T)
l1 = unique(spe$slide)
names(l1) = l1
l1 = lapply(l1, function(x) spe[,spe$slide==x])

for (i in seq_along(l1)) {
#	spe_small <- filter_genes(l1[[i]], filter_genes_ncounts = 3, filter_genes_pcspots = .5, filter_mito=T)
	spe_small = l1[[i]]
	cat("\n",names(l1)[i],"Calculating nnSVG... ",format(Sys.time(),tz="UTC"),"\n")
	cat("\n",dim(spe_small),"\n")
	spe_nnSVG <- nnSVG(spe_small, n_threads=12)
	svg = rowData(spe_nnSVG)
	cat("\nSaving results...\n\n")
	write.csv(svg, here("processed-data","04_preprocessing",paste0("nnSVG_",names(l1)[i],".csv")), row.names=T)
}

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
