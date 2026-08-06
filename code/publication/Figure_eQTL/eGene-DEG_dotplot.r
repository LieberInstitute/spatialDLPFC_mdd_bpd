setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
})
set.seed(123)

refined.modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv")

source("code/publication/plotting_utils.r")
de.df2 = filter(de.df, group!="MDD.BPD")

plot.genes = c("MAPK3","DUSP4","DUSP6","ELK1","TEF","RASD1",
	"FABP3","SURF1","ATP6V0E2",
	"SPON2","BAIAP3")

p3 <- getDotplot(plot.genes, de.df)#+theme(axis.text.y=element_text(size=7))
p3.1 <- getMeanRatioBar(plot.genes, sce_summ)

ggsave(file="plots/publication/Figure_eQTL/eQTL_dotplot.pdf", 
       arrangeGrob(grobs=list(p3, p3.1), layout_matrix=matrix(c(1,1,1,1,1,1,1,2), ncol=8)),
       height=3, width=6.5)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
