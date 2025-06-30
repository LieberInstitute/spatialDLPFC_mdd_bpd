setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
        library(Seurat)
	library(future)
})
set.seed(123)
cat("\nIncrease max size of global variables to fix error...\n")
options(future.globals.maxSize = 1000 * 1024^2) #1GB

cat("\nLoad files...\n")
load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT.Rdata")
load("processed-data/05_clustering/Seurat/seurat_MBv_processed-SCT.Rdata")

#remove low UMI spots
#cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
#stopifnot(identical(colnames(seu_mbv), rownames(cdata)))
#seu_mbv = seu_mbv[,cdata$precast_k9_1663_f!="low UMI"]
#cat("\nRemoved low UMI spots...\n")
#dim(seu_mbv)

#load in gene list
geneList = readRDS("processed-data/04_feature_selection/nnSVG-eval_geneList.rds")

##integrate
#cat("\nStart integration...\n")
#format(Sys.time())
#
#anchors <- FindTransferAnchors(reference = seu_con, 
#	features=intersect(VariableFeatures(seu_con),unlist(geneList)),
#	#features=geneList$qual_genes,
#	query = seu_mbv, normalization.method = "SCT",
#	k.anchor=50, k.score=50, max.features=500, npcs=20, dims=1:20, mapping.score.k=T)
#save(anchors, file="processed-data/05_clustering/Seurat/anchors_SZBDMulti-seq_MBv-filtered_ref-control_query-MBv_VarFeat-intersect-genes-kanchor-50-pc20.Rdata")
#cat("\nAnchors saved to: processed-data/05_clustering/Seurat/anchors_SZBDMulti-seq_MBv-filtered_ref-control_query-MBv_VarFeat-intersect-genes-kanchor-50-pc20.Rdata\n")

#load anchors
cat("\nLoad saved anchors...\n")
load("processed-data/05_clustering/Seurat/anchors_SZBDMulti-seq_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc20.Rdata")

cat("\nLoad in PRECAST results for query weights...\n")
format(Sys.time())
load("processed-data/05_clustering/PRECAST/srt_precast_k-9_n1663.Rdata")
#precast colnames are different and reordered
emb_mtx = seuInt[["PRECAST"]]@cell.embeddings
seurat_key = paste(seu_mbv$sample_id, colnames(seu_mbv), sep="_")
emb_mtx = emb_mtx[seurat_key,]
rownames(emb_mtx) = colnames(seu_mbv)
seu_mbv = PRECAST::Add_embed(emb_mtx, seu_mbv, embed_name="PRECAST", assay="SCT")

cat("\nTransfer labels...\n")
format(Sys.time())

mbv_query <- TransferData(anchorset = anchors, refdata = seu_con$seurat_annotated,
	#refdata = seurat_low.res,
	#weight.reduction = seu_mbv[["pca"]], dims = 1:30)
	weight.reduction = seu_mbv[["PRECAST"]], dims=1:15, k.weight=50)
write.csv(mbv_query, "processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc20_red-precast-kweight-50-seurat-annotated.csv", row.names=T)
cat("\nSaved label transfer to: processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc20_red-precast-kweight-50-seurat-annotated.csv\n")

## Reproducibility information
print("\n\nReproducibility information:")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
