setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
})
set.seed(123)

source("code/publication/plotting_utils.r")

p1 <- getViolin(c("EDN1","VEGFA"), version="seurat")
ggsave(file="plots/publication/Figure_ox-stress/EDN1-VEGFA_violin.pdf", p1, width=3.5, height=3.2)

p2 <- getViolin(c("APOLD1"), version="seurat")
ggsave(file="plots/publication/Figure_ox-stress/APOLD1_violin.pdf", p2, width=3.5, height=1.7)

p3 <- getViolin(c("A2M"), .y.upper.bound=12, .y.breaks=c(0,3,6,9,12), version="seurat")
ggsave(file="plots/publication/Figure_ox-stress/A2M_violin.pdf", p3, width=3.5, height=1.7)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
