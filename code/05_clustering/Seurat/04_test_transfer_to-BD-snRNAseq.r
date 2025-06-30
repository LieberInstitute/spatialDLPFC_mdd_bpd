setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(Seurat)
})
set.seed(123)

cat("\nLoad files...\n")
load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT.Rdata")
#if bipolar not yet processed
#load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_bipolar.Rdata")


##subset bipolar
##subset to only features that passed MBv pseudobulk gene filters
#cat("\nSubset to only genes that passed MBv pseudobulk gene filters...\n")
#var = read.csv("processed-data/05_clustering/Seurat/var_SZBDMulti-Seq_filtered.csv")
#stopifnot(identical(seu_bd[["RNA"]]@meta.data$featureid, var$featureid))
#seu_bd[["RNA"]]@meta.data$featurekey = var$featurekey
#cat("\nOriginal number of SZBDMulti-seq genes:",nrow(var),"\n")

#rdata = read.csv("processed-data/06_pseudobulk/pseudobulk-sample-n1663-k9_filtered-genes_avg-logcounts.csv", row.names=1)
#cat("\nOriginal number of MBv pseudobulk genes:",nrow(rdata),"\n")

#seu_bd = seu_bd[seu_bd[["RNA"]]@meta.data$featureid %in% rownames(rdata),]
#cat("\nNumber of MBv-filtered SZBDMulti-seq genes:",nrow(seu_bd),"\n")
#seu_bd

#process bipolar
#cat("\nPre-process bipolar dataset...\n")
#format(Sys.time())
#seu_bd <- SCTransform(seu_bd)
#seu_bd <- RunPCA(seu_bd)
#seu_bd <- FindNeighbors(seu_bd, dims = 1:20, k.param=50)
#seu_bd <- RunUMAP(seu_bd, dims=1:20)

#save(seu_bd, file="processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_bipolar_MBv-filtered_processed-SCT.Rdata")
#cat("\nBipolar processed and saved to: processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_bipolar_MBv-filtered_processed-SCT.Rdata\n")

#read in processed bipolar
load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_bipolar_MBv-filtered_processed-SCT.Rdata")

#load geneList
geneList = readRDS("processed-data/04_feature_selection/nnSVG-eval_geneList.rds")

##integrate
#cat("\nStart integration...\n")
#format(Sys.time())
#anchors <- FindTransferAnchors(reference = seu_con,
#       features=geneList$qual_genes,
#       query = seu_bd, normalization.method = "SCT",
#       k.anchor=50, k.score=50, max.features=500, npcs=20, dims=1:20, mapping.score.k=T)
#save(anchors, file="processed-data/05_clustering/Seurat/anchors_SZBDMulti-seq_MBv-filtered_ref-control_query-bipolar_qual-genes-kanchor-50-pc20.Rdata")
#cat("\nAnchors saved to: processed-data/05_clustering/Seurat/anchors_SZBDMulti-seq_MBv-filtered_ref-control_query-bipolar_qual-genes-kanchor-50-pc20.Rdata\n")

#read in anchors
cat("\nLoad saved anchors...\n")
load("processed-data/05_clustering/Seurat/anchors_SZBDMulti-seq_MBv-filtered_ref-control_query-bipolar_qual-genes-kanchor-50-pc20.Rdata")

cat("\nTransfer labels...\n")
format(Sys.time())

bd_query <- TransferData(anchorset = anchors, refdata = seu_con$azimuth,
        weight.reduction = seu_bd[["pca"]], dims = 1:30, k.weight=50)

write.csv(bd_query, "processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-bipolar_qual-genes-kanchor-50-pc20_red-pca-kweight-50-azimuth.csv", row.names=T)
cat("\nSaved label transfer to: processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-bipolar_qual-genes-kanchor-50-pc20_red-pca-kweight-50-azimuth.csv\n")


## Reproducibility information
print("Reproducibility information:")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
