library(SpatialExperiment)
library(ggspavis)
library(dplyr)
library(here)

load(here("processed-data","02_build_spe","spe_demo.Rdata"))

spe$remove = if_else(spe$in_tissue==FALSE, "off tissue","ok")
spe$remove = if_else(spe$sum_umi==0 & spe$in_tissue==TRUE, "zero counts", spe$remove)
spe$remove = factor(spe$remove, levels=c("ok","off tissue","zero counts"))

l1 = unique(spe$sample_id)
names(l1) = l1
l1 = lapply(l1, function(x) colData(spe)$sample_id==x)

p1 <- lapply(1:length(l1), function(x) plotSpots(spe[,l1[[x]]], annotate="remove")+scale_color_manual(values=c("black","grey","red"))+labs(title=names(l1)[[x]]))

pdf(here("plots", "03_QC", paste0("errant_spots_",Sys.Date(),".pdf")), width=8, height=12)
PRECAST::drawFigs(p1, layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
dev.off()

spe = spe[,spe$remove=="ok"]
save(spe, file=here("processed-data","03_QC","spe_demo-filt.Rdata"))
