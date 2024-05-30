setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(ggspavis)
	library(dplyr)
	library(here)
})

load(here("processed-data","02_build_spe","spe_demo.Rdata"))

spe$remove = if_else(spe$in_tissue==FALSE, "off tissue","ok")
spe$remove = if_else(spe$sum_umi==0 & spe$in_tissue==TRUE, "zero counts", spe$remove)
spe$remove = factor(spe$remove, levels=c("ok","off tissue","zero counts"))

l1 = unique(spe$sample_id)
names(l1) = lapply(l1, function(x) unique(colData(spe)[spe$sample_id==x,"brain"]))
l1 = lapply(l1, function(x) colData(spe)$sample_id==x)

p1 <- lapply(1:length(l1), function(x) plotSpots(spe[,l1[[x]]], annotate="remove",in_tissue=NULL)+scale_color_manual(values=c("black","darkgrey","red"))+labs(title=names(l1)[[x]]))

pdf(here("plots", "03_QC", paste0("errant_spots_",Sys.Date(),".pdf")), width=8, height=12)
PRECAST::drawFigs(p1, layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
dev.off()

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
