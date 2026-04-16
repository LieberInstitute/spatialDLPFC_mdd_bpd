suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(HDF5Array)
  library(dplyr)
  library(ggplot2)
  library(escheR)
  library(ggrastr)
  library(gridExtra)
})
set.seed(123)

cpList <- readRDS("plots/colorPalettes.rds")

#spot plots
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")
spotdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
spotdata = spotdata[rownames(colData(spe)),]
stopifnot(identical(rownames(colData(spe)), rownames(spotdata)))
spe$smoothed_k9_1663 = factor(spotdata$smoothed_k9_1663, levels=c("L1","L2","L3/4","L5","L6","WM","low UMI","Vasc","GABA"),
                              labels=c("L1","L2","L3.4","L5","L6","WM","dropped","dropped","dropped"))

res.pc30 <- read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc30_red-precast-kweight-50-low-res.csv", row.names=1)
res.pc30 = res.pc30[rownames(colData(spe)),]
stopifnot(identical(rownames(colData(spe)), rownames(res.pc30)))
spe$seurat_pc30 = factor(res.pc30$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","Inhb","L4","L5","L6","Oligo"),
                         labels=c("Micro.Vasc","Astro","L2.3","L2.3","Inhb","L4","L5","L6","Oligo"))


spe_sub = spe[,spe$sample_id=="V13B23-308_D1"]


spot.genes = c("FTL","APLP1","PLP1","ENPP2","WNK1","ANP32B","SGK1","NEAT1")
for(i in spot.genes) {
  spe_sub[[i]] = logcounts(spe_sub)[rowData(spe_sub)$gene_name==i,]
}

plist <- lapply(spot.genes, function(x) {
  p = make_escheR(spe_sub) %>%
    add_ground(var="smoothed_k9_1663", point_size=.5, stroke=.3) %>%
    add_fill(var=x, point_size=.6)
  p <- p+scale_color_manual(values=c(cpList$smoothed.light), guide="none")+
    scale_fill_gradient(low="white",high="black", guide="none")+
    labs(title=x)+
    theme(text=element_text(size=10), plot.title=element_text(face="italic"))
  return(rasterize(p, layers="Point", dpi=300))
})

plegend <- lapply(spot.genes, function(x) {
  p = make_escheR(spe_sub) %>%
    add_ground(var="smoothed_k9_1663", point_size=.1, stroke=.1) %>%
    add_fill(var=x, point_size=.1)
  p <- p+scale_color_manual(values=c(cpList$smoothed.light), guide="none")+
    scale_fill_gradient(low="white",high="black")+
    labs(title=x)+
    theme(text=element_text(size=10), plot.title=element_text(face="italic"))
  return(rasterize(p, layers="Point", dpi=300))
})

pdf(file="plots/publication/Figure3/supp_Oligo-mod-DEGs_spot-plot.pdf", width=4, height=8)
do.call(grid.arrange, c(plist, ncol=2))
do.call(grid.arrange, c(plegend, ncol=2))
dev.off()

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
