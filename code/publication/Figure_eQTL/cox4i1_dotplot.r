setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
})
set.seed(123)

refined.modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv")

source("code/publication/plotting_utils.r")

mod_gene="COX4I1"
plot.genes = c("ATP5F1E","NEDD8","COX7A1","SLC35E2B","TOMM7",
  "FABP3","PCP4","SURF1",
  "ATP6V0E2","NBL1",
  "ATF4","MAPK3")

p3 <- getDotplot(plot.genes, de.df)
p3.1 <- getMeanRatioBar(plot.genes, sce_summ)

ggsave(file=paste0("plots/publication/Figure_eQTL/",mod_gene,"_dotplot.pdf"), 
       arrangeGrob(grobs=list(p3, p3.1), layout_matrix=matrix(c(1,1,1,1,1,1,1,2), ncol=8)),
       height=3.2, width=6.5)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()

