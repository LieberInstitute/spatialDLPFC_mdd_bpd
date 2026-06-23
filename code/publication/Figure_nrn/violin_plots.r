setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
})
set.seed(123)

source("code/publication/plotting_utils.r")

#p1 <- getViolin(c("RASD1"), version="seurat")
#ggsave(file="plots/publication/Figure_nrn/RASD1_violin.pdf", p1, width=3.5, height=1.7)

#p1 <- getViolin(c("ELK1"), version="seurat")
#ggsave(file="plots/publication/Figure_nrn/ELK1_violin.pdf", p1, width=3.5, height=1.7)

#p2 <- getViolin(c("TEF"), version="seurat")
#ggsave(file="plots/publication/Figure_nrn/TEF_violin.pdf", p2, width=3.5, height=1.7)

p2 <- getViolin(c("TRBC2"), .y.upper.bound=9, .y.breaks=c(0,2,4,6,8), version="seurat")
p3 <- getViolin(c("TRBC2"), .y.upper.bound=9, .y.breaks=c(0,2,4,6,8), version="precast")
ggsave(file="plots/publication/Figure_nrn/TRBC2_violin.pdf", 
	marrangeGrob(grobs=list(p2, p3), ncol=1, nrow=1, top=NULL), 
	width=3.5, height=1.7)

#p3 <- getViolin(c("SST","CORT","CRH","VGF"), version="seurat")
#ggsave(file="plots/publication/Figure_nrn/SST-CORT-CRH-VGF_violin.pdf", p3, width=3.5, height=6)

p4 <- getViolin(c("PVALB","TAC1"), version="seurat")
p5 <- getViolin(c("PVALB","TAC1"), version="precast")
ggsave(file="plots/publication/Figure_nrn/supp_PVALB-TAC1_violin.pdf", 
	marrangeGrob(grobs=list(p4, p5), ncol=1, nrow=1, top=NULL),
	width=3.5, height=3.2)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
