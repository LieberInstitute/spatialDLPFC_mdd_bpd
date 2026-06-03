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



seu_con$azimuth_super.broad <- factor(as.character(seu_con$azimuth),
                                levels=c("Endo","PC","SMC","VLMC",
                                         "Immune","Micro",
                                         "Astro",
                                         "L2/3 IT", "L4 IT",
                                         "L5 IT", "L5 ET", "L5/6 NP",
                                         "L6 IT","L6 IT Car3","L6 CT","L6b",
                                         "Oligo","OPC",
                                         "Sncg","Pax6","Vip",
                                         "Lamp5","Lamp5 Lhx6",
                                         "Sst","Sst Chodl","Pvalb",
                                         "Chandelier"),
                                labels=c("Vasc","Vasc","Vasc","Vasc",
                                         "Micro","Micro",
                                         "Astro",
                                         "ExcN","ExcN",
                                         "ExcN","ExcN","ExcN",
                                         "ExcN","ExcN","ExcN","ExcN",
                                         "Oligo","Oligo",
                                         "InhN","InhN","InhN",
                                         "InhN","InhN",
                                         "InhN", "InhN", "InhN",
                                         "InhN")
)
color.palette = c("Vasc"=cpList$low.res.light[["Micro.Vasc"]],
            "Micro"=cpList$low.res.light[["L3"]],
            cpList$low.res.light[c("Astro")], cpList$low.res.bright["Oligo"],
            "InhN"=cpList$low.res.light[["Inhb"]],
            "ExcN"=cpList$low.res.light[["L2"]])

p2 <- DimPlot(seu_con, group.by="azimuth_super.broad",
        pt.size=.1, raster=F, label=F)+NoLegend()+
        scale_color_manual(values=color.palette)+
        theme(axis.title=element_blank(), axis.text=element_blank(), plot.title=element_blank())

ggsave(file="plots/publication/Figure1/seu-con_umap_super-broad-annotations.pdf",
        rasterize(p2, layers="Point", dpi=300),
        height=6, width=6)

load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_bipolar_MBv-filtered_processed-SCT.Rdata")
table(distinct(seu_bd@meta.data[,c("individualID","Biological_Sex")]))

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
