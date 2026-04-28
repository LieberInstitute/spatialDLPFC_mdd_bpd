setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(Seurat)
	library(dplyr)
	library(ggplot2)
	library(ggrastr)
})
set.seed(123)

load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT.Rdata")
seu_con

#head(seu_con@meta.data)
table(distinct(seu_con@meta.data[,c("individualID","Biological_Sex")]))

cpList <- readRDS("plots/colorPalettes.rds")
seu_con@meta.data$seurat_low.res = factor(as.character(seu_con@meta.data$seurat_low.res), 
                                          levels=c("Micro/Vasc","Astro","Oligo","L2","L3","L4","L5","L6","Inhb"),
                                          labels=c("Micro.Vasc","Astro","Oligo","L2","L3","L4","L5","L6","Inhb"))

p1 <- DimPlot(seu_con, group.by="seurat_low.res",
	pt.size=.1, raster=F, label=F)+NoLegend()+
	scale_color_manual(values=cpList$low.res.bright)+
	theme(axis.title=element_blank(), axis.text=element_blank(), plot.title=element_blank())

ggsave(file="plots/publication/Figure1/seu-con_umap_cell-type-annotations.pdf",
	rasterize(p1, layers="Point", dpi=300),
	height=6, width=6)



load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_bipolar_MBv-filtered_processed-SCT.Rdata")
table(distinct(seu_bd@meta.data[,c("individualID","Biological_Sex")]))

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
