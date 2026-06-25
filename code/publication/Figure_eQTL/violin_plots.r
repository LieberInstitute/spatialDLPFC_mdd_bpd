setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
})
set.seed(123)

source("code/publication/plotting_utils.r")

# these are all part of COX4I1 module
#p1 <- getViolin(c("SURF1","ATP6V0E2","FABP3"),version="seurat")
#ggsave(file="plots/publication/Figure_eQTL/SURF1-ATP6V0E2-FABP3_violin.pdf", p1, width=3.5, height=4.3)

#p1 <- getViolin(c("STMN4"), .y.upper.bound=13, .y.breaks=c(0,3,6,9,12), version="seurat")
#ggsave(file="plots/publication/Figure_eQTL/STMN4_violin.pdf", p1, width=3.5, height=1.7)

#p2 <- getViolin(c("IFITM3","IFITM2"), .y.upper.bound=13, .y.breaks=c(0,3,6,9,12), version="seurat")
#ggsave(file="plots/publication/Figure_eQTL/IFITM3-IFITM2_violin.pdf", p2, width=3.5, height=3.2)

#p3 <- getViolin(c("MAPK3","DUSP6","DUSP4"), version="seurat")
#ggsave(file="plots/publication/Figure_eQTL/MAPK3-DUSP6-DUSP4_violin.pdf", p3, width=3.5, height=4.3)

p3 <- getViolin(c("DUSP6"), version="seurat")
ggsave(file="plots/publication/Figure_eQTL/DUSP6_violin.pdf", p3, width=3.5, height=1.7)

plist4 <- getViolin(c("MAPK3","DUSP4"), .y.upper.bound=10, .y.breaks=c(0,3,6,9), version="seurat", whole.tissue_only=T)
ggsave(file="plots/publication/Figure_eQTL/MAPK3-DUSP4_whole-tissue-only_violin.pdf", 
	arrangeGrob(grobs=plist4, top=NULL, ncol=2), width=2, height=1.7)


## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
