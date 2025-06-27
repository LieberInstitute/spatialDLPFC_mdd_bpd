setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(Seurat)
})
set.seed(123)

load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control.Rdata")

#subset to only features that passed MBv pseudobulk gene filters
cat("\nSubset to only genes that passed MBv pseudobulk gene filters...\n")
var = read.csv("processed-data/05_clustering/Seurat/var_SZBDMulti-Seq_filtered.csv")
stopifnot(identical(seu_con[["RNA"]]@meta.data$featureid, var$featureid))
seu_con[["RNA"]]@meta.data$featurekey = var$featurekey
cat("\nOriginal number of SZBDMulti-seq genes:",nrow(var),"\n")

rdata = read.csv("processed-data/06_pseudobulk/pseudobulk-sample-n1663-k9_filtered-genes_avg-logcounts.csv", row.names=1)
cat("\nOriginal number of MBv pseudobulk genes:",nrow(rdata),"\n")

seu_con = seu_con[seu_con[["RNA"]]@meta.data$featureid %in% rownames(rdata),]
cat("\nNumber of MBv-filtered SZBDMulti-seq genes:",nrow(seu_con),"\n")
seu_con


cat("\n\nSCT and processing...\n")
format(Sys.time())
seu_con <- SCTransform(seu_con, ncells=5000)
#seu_con <- NormalizeData(seu_con)
#seu_con <- FindVariableFeatures(seu_con)
#seu_con <- ScaleData(seu_con)

save(seu_con, file="processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT.Rdata")
cat("\nSaved checkpoint file to: processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT.Rdata\n\n")

seu_con <- RunPCA(seu_con)
seu_con <- RunUMAP(seu_con, dims=1:20)

seu_con <- FindNeighbors(seu_con, dims = 1:20, k.param=50) 
seu_con <- FindClusters(seu_con, resolution=.2) 

seu_con

save(seu_con, file="processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT.Rdata")
cat("\nSaved to: processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT.Rdata\n")

print("Reproducibility information:")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
