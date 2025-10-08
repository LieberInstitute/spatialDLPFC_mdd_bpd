setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(Seurat)
	library(SpatialExperiment)
	library(ggplot2)
})
set.seed(123)

#load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control.Rdata")

##subset to only features that passed MBv pseudobulk gene filters
#cat("\nSubset to only genes that passed MBv pseudobulk gene filters...\n")
#var = read.csv("processed-data/05_clustering/Seurat/var_SZBDMulti-Seq_filtered.csv")
#stopifnot(identical(seu_con[["RNA"]]@meta.data$featureid, var$featureid))
#seu_con[["RNA"]]@meta.data$featurekey = var$featurekey
#cat("\nOriginal number of SZBDMulti-seq genes:",nrow(var),"\n")

#load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_conservative_pseudo_sample-smoothed-n1626-k7_norm-filt.Rdata")
#pass.genes = rowSums(as.data.frame(rowData(spe_pseudo)[,c("high_expr_group_sample_id","high_expr_group_cluster","high_expr_SZBDMultiseq")]))>2

##rdata = read.csv("processed-data/06_pseudobulk/PRECAST_smoothed/pseudobulk-sample-smoothed-n1626-k7_conservative_filtered-genes_avg-logcounts.csv", row.names=1)
##cat("\nOriginal number of MBv pseudobulk genes:",nrow(rdata),"\n")

#seu_con = seu_con[seu_con[["RNA"]]@meta.data$featureid %in% rownames(spe_pseudo)[pass.genes],]
#cat("\nNumber of MBv-filtered SZBDMulti-seq genes:",nrow(seu_con),"\n")
#seu_con


#cat("\n\nSCT and processing...\n")
#format(Sys.time())
#seu_con <- SCTransform(seu_con, ncells=5000)

#save(seu_con, file="processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered-conservative_processed-SCT.Rdata")
#cat("\nSaved checkpoint file to: processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered-conservative_processed-SCT.Rdata\n\n")
cat("\nLoad saved checkpoint file: processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered-conservative_processed-SCT.Rdata\n")
load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered-conservative_processed-SCT.Rdata")
seu_con

seu_con <- RunPCA(seu_con)
seu_con <- RunUMAP(seu_con, dims=1:15)

seu_con <- FindNeighbors(seu_con, dims = 1:20, k.param=50) 
#load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT.Rdata")
seu_con <- FindClusters(seu_con, resolution=.3) 

# additional cell type annotations
cat("\nAzimuth-based annotations...\n")
seu_con$azimuth_broad <- factor(as.character(seu_con$azimuth),
        levels=c("Endo","Immune","Micro","PC","SMC","VLMC",
                "Astro",
                "L2/3 IT", "L4 IT",
                "L5 IT", "L5 ET", "L5/6 NP",
                "L6 IT","L6 IT Car3","L6 CT","L6b",
                "Oligo","OPC",
                "Chandelier","Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Sncg","Sst","Sst Chodl","Vip"),
        labels=c("Micro/Vasc","Micro/Vasc","Micro/Vasc","Micro/Vasc","Micro/Vasc","Micro/Vasc",
                "Astro",
                "L2/3","L4",
                "L5","L5","L5",
                "L6","L6","L6","L6",
                "Oligo","Oligo",
                "Inhb","Inhb","Inhb","Inhb","Inhb","Inhb","Inhb","Inhb","Inhb")
)
table(cbind.data.frame("azimuth"=seu_con$azimuth, "azimuth_broad"=seu_con$azimuth_broad))

cat("\nSeurat processing-based annotations...\n")
table(seu_con@meta.data[,c("azimuth","seurat_clusters")])

#azmap = c("L3"="0","L2"="1","L5"="2","Oligo"="3","L4/L5"="4","L6"="5",
#          "Inhb CGE VIP"="6","Inhb MGE PV"="7", "L6b"="8","Inhb MGE SST"="9","Astro"="10",
#          "L3/L4"="11","Inhb CGE LAMP5"="12","OPC"="13","L5/6 NP"="14","Micro/Vasc"="15",
#          "L6 IT Car3"="16", "Chandelier"="17", "L5 ET"="18")
#seu_con$seurat_annotated = factor(as.character(seu_con$seurat_clusters), levels=azmap,
#                                  labels=names(azmap))

#seu_con$seurat_low.res <- factor(as.character(seu_con$seurat_annotated),
#                                 levels=c("Micro/Vasc","Astro",
#                                          "Inhb MGE PV","Inhb MGE SST","Inhb CGE VIP","Inhb CGE LAMP5","Chandelier",
#                                          "L2","L3",
#                                          "L3/L4","L4/L5",
#                                          "L5","L5/6 NP","L5 ET",
#                                          "L6","L6 IT Car3","L6b",
#                                          "Oligo","OPC"),
#                                 labels=c("Micro/Vasc","Astro",
#                                          "Inhb","Inhb","Inhb","Inhb","Inhb",
#                                          "L2","L3",
#                                          "L4","L4",
#                                          "L5","L5","L5",
#                                          "L6","L6","L6",
#                                          "Oligo","Oligo"))
#table(cbind.data.frame("annotated_clusters"=seu_con$seurat_annotated, "low.res_clusters"=seu_con$seurat_low.res))

azmap = c("L3"=0, "L5"=1, "L2"=2, "L4/L5"=3, "L6"=4, "Inhb CGE VIP"=5, "Oligo"=6,
	"Inhb MGE PV"=7, "L6b"=8, "Inhb MGE SST"=9, "Oligo"=10, "Astro"=11, "L3/L4"=12,
	"Inhb CGE LAMP5"=13, "OPC"=14, "Micro/Vasc"=15, "L5/6 NP"=16, "L6 IT Car3"=17, 
	"Chandelier"=18, "L5 ET"=19)
seu_con$seurat_annotated = factor(as.character(seu_con$seurat_clusters), levels=azmap,
                                  labels=names(azmap))
seu_con$seurat_low.res <- factor(as.character(seu_con$seurat_annotated),
                                 levels=c("Micro/Vasc","Astro",
                                          "Inhb MGE PV","Inhb MGE SST","Inhb CGE VIP","Inhb CGE LAMP5","Chandelier",
                                          "L2","L3",
                                          "L3/L4","L4/L5",
                                          "L5","L5/6 NP","L5 ET",
                                          "L6","L6 IT Car3","L6b",
                                          "Oligo","Oligo","OPC"),
                                 labels=c("Micro/Vasc","Astro",
                                          "Inhb","Inhb","Inhb","Inhb","Inhb",
                                          "L2","L3",
                                          "L4","L4",
                                          "L5","L5","L5",
                                          "L6","L6","L6",
                                          "Oligo","Oligo","Oligo"))
table(cbind.data.frame("annotated_clusters"=seu_con$seurat_annotated, "low.res_clusters"=seu_con$seurat_low.res))

save(seu_con, file="processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered-conservative_processed-SCT.Rdata")
cat("\nSaved to: processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered-conservative_processed-SCT.Rdata\n")

#plot annotations
cat("\nUMAP plots of annotations...\n")
cpList <- readRDS("plots/colorPalettes.rds")
seu_con@meta.data$seurat_low.res = factor(as.character(seu_con@meta.data$seurat_low.res), 
                                          levels=c("Micro/Vasc","Astro","Oligo","L2","L3","L4","L5","L6","Inhb"),
                                          labels=c("Micro.Vasc","Astro","Oligo","L2","L3","L4","L5","L6","Inhb"))
p1 <- DimPlot(seu_con, group.by=c("azimuth","azimuth_broad",#"seurat_clusters"),
	"seurat_annotated","seurat_low.res"),
	pt.size=.1, raster=F, label=T, combine=F)
ggsave(file="plots/05_clustering/Seurat/seu-con_conservative_umap_cell-type-annotations.png",
	gridExtra::grid.arrange(p1[[1]]+NoLegend(), p1[[2]]+NoLegend(),
		p1[[3]]+NoLegend(), 
		p1[[4]]+scale_color_manual(values=cpList$low.res.bright)+NoLegend(),
		ncol=2),
	bg="white", height=16, width=12)
cat("\nSaved to: plots/05_clustering/Seurat/seu-con_conservative_umap_cell-type-annotations.png\n")

print("Reproducibility information:")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
