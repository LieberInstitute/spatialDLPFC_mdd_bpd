setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(scater)
  library(gridExtra)
  library(ggplot2)
  library(gridExtra)
  library(ggrastr)
})

set.seed(123)


cpList = readRDS("plots/colorPalettes.rds")

#PRECAST smoothed
load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")

p1 <- plotReducedDim(spe_pseudo, dimred="PCA_1663", ncomponents=2, colour_by = "smoothed_k9_1663", point_alpha=1)+
  scale_color_manual("", values=cpList$smoothed.bright)+
p2 <- plotReducedDim(spe_pseudo, dimred="PCA_1663", ncomponents=2, colour_by = "condition", shape_by= "sex", point_alpha=1)+
  scale_color_manual("", values=cpList$dx.pal)+
  scale_shape_manual(values=c(21,23))
b1 = ggplot_build(p1)
b1$data[[1]]$size = .5
b1$data[[1]]$shape = 21
# because of overlap of points I like empty fill better
b1$data[[1]]$fill = NA

b2 = ggplot_build(p2)
b2$data[[1]]$size = .5
b2$data[[1]]$fill = NA
## alt to fill with lighter version of dx.pal color
## based on the small size of the points I actually think i like empty better
#b2$data[[1]]$fill = as.character(factor(b2$data[[1]]$colour, levels=c(cpList$dx.pal[[1]], cpList$dx.pal[[2]], cpList$dx.pal[[3]]),
#                             labels=c(colorRampPalette(c("white",cpList$dx.pal[[1]]))(20)[[12]],
#                                      colorRampPalette(c("white",cpList$dx.pal[[2]]))(20)[[12]],
#                                      colorRampPalette(c("white",cpList$dx.pal[[3]]))(20)[[12]])))

ggsave(file="plots/publication/Figure1/PC1-PC2_UMAP.pdf",
       grid.arrange(ggplot_gtable(b1), ggplot_gtable(b2), ncol=1),
       width=3, height=4)

