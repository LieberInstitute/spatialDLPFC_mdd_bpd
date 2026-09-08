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

#load precast clusters
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
stopifnot(identical(rownames(colData(spe)), rownames(cdata)))

spe$smoothed = factor(cdata$smoothed_k9_1663_f, levels=c("L1","L2","L3/4","L5","L6","WM","Vasc","GABA","low UMI"),
                     labels=c("L1","L2","L3.4","L5","L6","WM","Vasc","GABA","low UMI"))


#load seurat clusters
res = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc30_red-precast-kweight-50-low-res.csv", row.names=1)
stopifnot(identical(rownames(colData(spe)), rownames(res)))

spe$merged = factor(res$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"),
	labels=c("Micro.Vasc","Astro","L2.3","L2.3","L4","Inhb","L5","L6","Oligo"))


spe_sub = spe[,spe$sample_id=="V13B23-308_D1"]

table(colData(spe_sub)[,c("smoothed","merged")])

# domain-SP spot plot
spe_sub2 = spe_sub[,!spe_sub$smoothed %in% c("Vasc","GABA","low UMI")]
p = make_escheR(spe_sub2) %>%
  add_fill(var="smoothed", point_size = .9)
p1 <- p+scale_fill_manual(values=cpList$smoothed.bright, guide="none")+
  theme(plot.title=element_text(size=14, hjust=.5))

# domain-CT spot plot
p = make_escheR(spe_sub) %>%
  add_fill(var="merged", point_size = .9)
p2 <- p+scale_fill_manual(values=cpList$transfer.bright, guide="none")+
  theme(plot.title=element_text(size=14, hjust=.5))



pdf(file="plots/publication/Figure1/domain-SP-CT_spot-plots.pdf", width=4, height=2.5)
grid.arrange(rasterize(p1, dpi=300), rasterize(p2, dpi=300), ncol=2)
dev.off()

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
