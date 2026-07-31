setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
})
set.seed(123)

source("code/publication/plotting_utils.r")

#print(filter(sig.df, gene_name=="ANP32B", source=="se", sex.group=="MDD.BPD"))
#p1 <- getViolin(c("ANP32B"), version="seurat")
p1 <- getViolin(c("ANP32B"), version="precast")

ggsave(file="plots/publication/Figure3/supp_ANP32B_violin-precast.pdf", p1, width=3.5, height=1.7)

#print(filter(sig.df, gene_name %in% c("CLDN11","SGK1"), source=="se", sex.group=="MDD.BPD"))
#p2 <- getViolin(c("CLDN11","SGK1"), .y.upper.bound=12, .y.breaks=c(0,3,6,9,12), version="seurat")
p2 <- getViolin(c("CLDN11","SGK1"), .y.upper.bound=12, .y.breaks=c(0,3,6,9,12), version="precast")
ggsave(file="plots/publication/Figure3/supp_CLDN11-SGK1_violin-precast.pdf", p2, width=3.5, height=3.2)


#print(filter(sig.df, gene_name %in% c("BCL6","HSPA1B","TPT1","UBA52"), source=="se", sex.group=="MDD.BPD"))
#p3 <- getViolin(c("BCL6","HSPA1B","TPT1","UBA52"), .y.upper.bound=14, .y.breaks=c(0,3,6,9,12), version="seurat")
#ggsave(file="plots/publication/Figure3/supp_BCL6-HSPA1B-TPT1-UBA52_violin.pdf", p3, width=3.5, height=7)


## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
