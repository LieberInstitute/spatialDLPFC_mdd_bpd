setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
})
set.seed(123)

source("code/publication/plotting_utils.r")

p1 <- getViolin(c("EDN1","VEGFA"), version="seurat")
ggsave(file="plots/publication/Figure4/EDN1-VEGFA_violin.pdf", p1, width=3.5, height=3.2)

p2 <- getViolin(c("APOLD1"), version="seurat")
ggsave(file="plots/publication/Figure4/APOLD1_violin.pdf", p2, width=3.5, height=1.7)


#p4 <- getViolin(c("APOLD1","ANGPTL4"), #.y.upper.bound=13, .y.breaks=c(0,3,6,9,12), 
#	version="seurat")
#ggsave(file="plots/publication/Figure4/APOLD1-ANGPTL4_violin.pdf", p4, width=3.5, height=3.2)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
