setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(scater)
  library(gridExtra)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(ggrastr)
})

set.seed(123)


cpList = readRDS("plots/colorPalettes.rds")

load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")

p1 <- plotReducedDim(spe_pseudo, dimred="PCA_1663", ncomponents=2, colour_by = "seurat_label", point_alpha=1)+
  scale_color_manual("", values=cpList$transfer.bright)
p2 <- plotReducedDim(spe_pseudo, dimred="PCA_1663", ncomponents=2, colour_by = "condition", shape_by= "sex", point_alpha=1)+
  scale_color_manual("", values=cpList$dx.pal)+
  scale_shape_manual(values=c(21,23))
b1 = ggplot_build(p1+theme(legend.position="none"))
b1$data[[1]]$size = .5
b1$data[[1]]$shape = 21
# because of overlap of points I like empty fill better
b1$data[[1]]$fill = NA

b2 = ggplot_build(p2+theme(legend.position="none"))
b2$data[[1]]$size = .5
b2$data[[1]]$fill = NA

pdf("plots/publication/supp_clustering_Seurat/PC1-PC2.pdf", width= 2, height=2)
plot(ggplot_gtable(b1))
p1
plot(ggplot_gtable(b2))
p2
dev.off()

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
