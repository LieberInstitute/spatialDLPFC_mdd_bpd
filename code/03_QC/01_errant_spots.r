setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(dplyr)
	library(here)
})

load(here("processed-data","02_build_spe","spe_demo.Rdata"))

spe$remove = if_else(spe$in_tissue==FALSE, "off tissue","ok")
spe$remove = if_else(spe$sum_umi==0 & spe$in_tissue==TRUE, "zero counts", spe$remove)
spe$remove = factor(spe$remove, levels=c("ok","off tissue","zero counts"))

spe = spe[,spe$remove=="ok"]
save(spe, file=here("processed-data","03_QC","spe_demo-filt.Rdata"))

write(c(paste("Created spe_demo-filt on",format(Sys.time(), tz="UTC"),"UTC"),
        paste("File location:",here("processed-data","03_QC","spe_demo-filt.Rdata")),
        paste("Source code:",here("code","03_QC","01_errant_spots.r")),
        "*","*","*"), here("spe_tracker_current.txt"), append=TRUE)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
