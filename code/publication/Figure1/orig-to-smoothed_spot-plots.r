setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(HDF5Array)
  library(DelayedArray)
  library(dplyr)
  library(ggplot2)
  library(escheR)
  library(gridExtra)
  library(ggrastr)
})


set.seed(123)
setAutoBlockSize(1e9)

cpList = readRDS("plots/colorPalettes.rds")

#load spe
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")

#load clusters
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
stopifnot(identical(rownames(colData(spe)), rownames(cdata)))

colnames(cdata)

spe$precast = factor(cdata$precast_k9_1663_f, levels=c("L1","L2","L3/4","L5","L6","WM","Vasc","GABA","low UMI"),
                     labels=c("L1","L2","L3.4","L5","L6","WM","Vasc","GABA","low UMI"))
spe$smoothed = factor(cdata$smoothed_k9_1663_f, levels=c("L1","L2","L3/4","L5","L6","WM","Vasc","GABA","low UMI"),
                     labels=c("L1","L2","L3.4","L5","L6","WM","Vasc","GABA","low UMI"))


spe_sub = spe[,spe$sample_id=="V13B23-308_D1"]

table(colData(spe_sub)[,c("precast","smoothed")])

names(cpList)
names(cpList$earthy.pal2)


p = make_escheR(spe_sub) %>%
  add_fill(var="precast", point_size = .8)
p1 <- p+scale_fill_manual(values=c(cpList$earthy.pal2, "low UMI"="grey"), guide="none")+
  theme(plot.title=element_text(size=14, hjust=.5))

p = make_escheR(spe_sub) %>%
  add_fill(var="smoothed", point_size = .8)
p2 <- p+scale_fill_manual(values=c(cpList$earthy.pal2, "low UMI"="grey"), guide="none")+
  theme(plot.title=element_text(size=14, hjust=.5))

spe_sub2 = spe_sub[,!spe_sub$smoothed %in% c("Vasc","GABA","low UMI")]
p = make_escheR(spe_sub2) %>%
  add_fill(var="smoothed", point_size = .8)
p3 <- p+scale_fill_manual(values=cpList$earthy.pal2, guide="none")+
  theme(plot.title=element_text(size=14, hjust=.5))

pdf(file="plots/publication/Figure1/supp_original-smoothed_spot-plot.pdf", width=6, height=2.5)
grid.arrange(rasterize(p1, dpi=300), rasterize(p2, dpi=300), rasterize(p3, dpi=300), ncol=3)
dev.off()


## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
