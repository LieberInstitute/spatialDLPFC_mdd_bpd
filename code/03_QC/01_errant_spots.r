library(SpatialExperiment)
library(dplyr)
library(here)

load(here("processed-data","02_build_spe","spe_demo.Rdata"))

spe$remove = if_else(spe$in_tissue==FALSE, "off tissue","ok")
spe$remove = if_else(spe$sum_umi==0 & spe$in_tissue==TRUE, "zero counts", spe$remove)
spe$remove = factor(spe$remove, levels=c("ok","off tissue","zero counts"))

spe = spe[,spe$remove=="ok"]
save(spe, file=here("processed-data","03_QC","spe_demo-filt.Rdata"))
