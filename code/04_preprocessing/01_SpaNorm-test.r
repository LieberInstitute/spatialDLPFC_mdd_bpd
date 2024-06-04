setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(SpaNorm)
	library(here)
})
set.seed(123)

load(here("processed-data","03_QC","spe_demo-filt.Rdata"))
spe = nnSVG::filter_genes(spe, filter_genes_ncounts = 3, filter_genes_pcspots = .5)
#single sample first just to see how long it takes
#spe_small = spe[,spe$brain=="Br6529"]
#then single slide
#spe_small = spe[,spe$brain %in% unique(spe$brain)[c(1,5,9,13,17,21)]]
#then whole thing
spe = SpaNorm(spe, verbose=TRUE)
save(spe, file=here("processed-data","04_preprocessing","spe_SpaNorm.Rdata"))

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
