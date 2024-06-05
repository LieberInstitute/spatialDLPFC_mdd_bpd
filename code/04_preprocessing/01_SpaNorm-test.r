setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(SpaNorm)
	library(here)
})
set.seed(123)

load(here("processed-data","03_QC","spe_demo-filt.Rdata"))

# previously used nnSVG filter to limit genes, now using feature list
#spe = nnSVG::filter_genes(spe, filter_genes_ncounts = 3, filter_genes_pcspots = .5)
feature.list = readRDS(here("processed-data","04_preprocessing","bindev-2k-3k_svg_feature-list.rda"))
#single sample first just to see how long it takes
#spe_small = spe[,spe$brain=="Br6529"]
#then 6 samples
spe_small = spe[feature.list$bindev.3k,spe$brain %in% unique(spe$brain)[c(1,5,9,13,17,21)]]
#then whole thing

spe_small = SpaNorm(spe_small, scale.factor=2, verbose=TRUE)
save(spe_small, file=here("processed-data","04_preprocessing","spe_SpaNorm-sf2-test6-bindev.3k.Rdata"))

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
