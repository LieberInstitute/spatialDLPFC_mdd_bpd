setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
})
set.seed(123)

refined.modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv")

source("code/publication/plotting_utils.r")

calcHeight <- function(ngenes) round(0.08*ngenes+2.7, 1)

for(i in unique(refined.modules$TF)) {
	plot.genes = c(i, filter(refined.modules, TF==i)$target)
	
	p3 <- getDotplot(plot.genes, de.df)
	p3.1 <- getMeanRatioBar(plot.genes, sce_summ)

	ggsave(paste0("plots/publication/modules/all_components/", i, "_dotplot.pdf"),
		arrangeGrob(grobs=list(p3, p3.1), layout_matrix=matrix(c(1,1,1,1,1,1,1,2), ncol=8), top=i),
		height=calcHeight(length(plot.genes)), width=6.5)
}

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
