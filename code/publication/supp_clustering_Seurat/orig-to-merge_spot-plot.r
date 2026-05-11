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
res = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc30_red-precast-kweight-50-low-res.csv", row.names=1)
stopifnot(identical(rownames(colData(spe)), rownames(res)))

spe$original = factor(res$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"),
	labels=c("Micro.Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"))
spe$merged = factor(res$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"),
	labels=c("Micro.Vasc","Astro","L2.3","L2.3","L4","Inhb","L5","L6","Oligo"))


spe_sub = spe[,spe$sample_id=="V13B23-308_D1"]

table(colData(spe_sub)[,c("original","merged")])


p = make_escheR(spe_sub) %>%
  add_fill(var="original", point_size = .8)
p1 <- p+scale_fill_manual(values=cpList$low.res.bright, guide="none")+
  theme(plot.title=element_text(size=14, hjust=.5))

p = make_escheR(spe_sub) %>%
  add_fill(var="merged", point_size = .8)
p2 <- p+scale_fill_manual(values=cpList$transfer.bright, guide="none")+
  theme(plot.title=element_text(size=14, hjust=.5))

pdf(file="plots/publication/supp_clustering_Seurat/original-merged_spot-plot.pdf", width=6, height=2.5)
grid.arrange(rasterize(p1, dpi=300), rasterize(p2, dpi=300), rasterize(p2, dpi=300), ncol=3)
dev.off()


## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
