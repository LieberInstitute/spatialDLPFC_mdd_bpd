setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
        library(Seurat)
	library(future)
})
set.seed(123)
cat("\nIncrease max size of global variables to fix error...\n")
options(future.globals.maxSize = 1000 * 1024^2) #1GB

cat("\nLoad files...\n")
load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered-conservative_processed-SCT.Rdata")
cat("\n\nSZBDMultiseq control:\n\n")
seu_con

load("processed-data/05_clustering/Seurat/seurat_MBv_conservative_processed-SCT.Rdata")
cat("\n\nMBv conservative:\n\n")
seu_mbv

##remove low UMI spots
##cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
##stopifnot(identical(colnames(seu_mbv), rownames(cdata)))
##seu_mbv = seu_mbv[,cdata$precast_k9_1663_f!="low UMI"]
##cat("\nRemoved low UMI spots...\n")
##dim(seu_mbv)

#subset to MBv spots with precast embeddings before finding anchors
cat("\nLoad in PRECAST results for query weights...\n")
format(Sys.time())
load("processed-data/05_clustering/PRECAST/srt_precast-conservative_k-7_n1626.Rdata")
#precast colnames are different and reordered
emb_mtx = seuInt[["PRECAST"]]@cell.embeddings
cat("\nDim PRECAST embedding matrix:\n")
dim(emb_mtx)

if(dim(emb_mtx)[1]!=dim(seu_mbv)[2]) {
        cat("\nDims of PRECAST embeddings and seu_mbv not equal\n")
        #seu_mbv$seurat_key = paste(seu_mbv$sample_id, colnames(seu_mbv), sep="_")
        both.spots = intersect(seu_mbv$seurat_key, rownames(emb_mtx))
        cat("\nNumber of spots with PRECAST embeddings:", length(both.spots), "\n")
        emb_mtx = emb_mtx[both.spots,]
        #rownames(emb_mtx) = both.spots
        seu_mbv = seu_mbv[,seu_mbv$seurat_key %in% both.spots]
        #cat("\nDim filtered seu_mbv to all spots with PRECAST embeddings:", dim(seu_mbv), "\n")
        #emb_mtx = emb_mtx[seu_mbv$seurat_key,]
        ##stopifnot(identical(rownames(emb_mtx), seu_mbv$seurat_key))
        #### for some reason identical is not working but they are all the same??
        #stopifnot(table(rownames(emb_mtx)==seu_mbv$seurat_key)[["TRUE"]]==length(both.spots))
        #rownames(emb_mtx) = colnames(seu_mbv)
        #seu_mbv = PRECAST::Add_embed(emb_mtx, seu_mbv, embed_name="PRECAST", assay="SCT")
}

#load in gene list
geneList = readRDS("processed-data/04_feature_selection/nnSVG-eval_conservative_geneList.rds")

#integrate
cat("\nStart integration...\n")
format(Sys.time())

anchors <- FindTransferAnchors(reference = seu_con, 
	#features=intersect(VariableFeatures(seu_con),unlist(geneList)),
	features=geneList$qual_genes,
	query = seu_mbv, normalization.method = "SCT",
	k.anchor=50, k.score=50, max.features=500, npcs=20, dims=1:20, mapping.score.k=T)
save(anchors, file="processed-data/05_clustering/Seurat/anchors_SZBDMulti-seq_MBv-filtered-conservative_ref-control_query-MBv-conservative_qual-genes-kanchor-50-pc20.Rdata")
cat("\nAnchors saved to: processed-data/05_clustering/Seurat/anchors_SZBDMulti-seq_MBv-filtered-conservative_ref-control_query-MBv-conservative_qual-genes-kanchor-50-pc20.Rdata\n")

##load anchors
#cat("\n\nLoad saved anchors...\n")
#load("processed-data/05_clustering/Seurat/anchors_SZBDMulti-seq_MBv-filtered-conservative_ref-control_query-MBv-conservative_qual-genes-kanchor-50-pc30.Rdata")

cat("\nLoad in PRECAST results for query weights...\n")
format(Sys.time())
load("processed-data/05_clustering/PRECAST/srt_precast-conservative_k-7_n1626.Rdata")
#precast colnames are different and reordered
emb_mtx = seuInt[["PRECAST"]]@cell.embeddings
cat("\nDim PRECAST embedding matrix:\n")
dim(emb_mtx)

#seurat_key = paste(seu_mbv$sample_id, colnames(seu_mbv), sep="_")
emb_mtx = emb_mtx[seu_mbv$seurat_key,]
stopifnot(table(rownames(emb_mtx)==seu_mbv$seurat_key)[["TRUE"]]==dim(emb_mtx)[1])
rownames(emb_mtx) = colnames(seu_mbv)
seu_mbv = PRECAST::Add_embed(emb_mtx, seu_mbv, embed_name="PRECAST", assay="SCT")

#cat("\nDim seu_mbv:\n")
#dim(seu_mbv)

#stopifnot(identical(dim(emb_mtx)[1], dim(seu_mbv)[2]))
#if(dim(emb_mtx)[1]!=dim(seu_mbv)[2]) {
#	cat("\nDims of PRECAST embeddings and seu_mbv not equal\n")
#	#seu_mbv$seurat_key = paste(seu_mbv$sample_id, colnames(seu_mbv), sep="_")
#	both.spots = intersect(seu_mbv$seurat_key, rownames(emb_mtx))
#	cat("\nNumber of spots with PRECAST embeddings:", length(both.spots), "\n")
#	emb_mtx = emb_mtx[both.spots,]
#	#rownames(emb_mtx) = both.spots
#	seu_mbv = seu_mbv[,seu_mbv$seurat_key %in% both.spots]
#	cat("\nDim filtered seu_mbv to all spots with PRECAST embeddings:", dim(seu_mbv), "\n")
#	emb_mtx = emb_mtx[seu_mbv$seurat_key,]
#	#stopifnot(identical(rownames(emb_mtx), seu_mbv$seurat_key))
#	### for some reason identical is not working but they are all the same??
#	stopifnot(table(rownames(emb_mtx)==seu_mbv$seurat_key)[["TRUE"]]==length(both.spots))
#	rownames(emb_mtx) = colnames(seu_mbv)
#	seu_mbv = PRECAST::Add_embed(emb_mtx, seu_mbv, embed_name="PRECAST", assay="SCT")
#}
cat("\nseu_mbv with PRECAST embeddings:\n\n")
seu_mbv

cat("\nCheck labels...\n")
table(seu_con$seurat_low.res, useNA="ifany")

cat("\n\nTransfer labels...\n")
format(Sys.time())

mbv_query <- TransferData(anchorset = anchors, #refdata = seu_con$seurat_annotated,
	refdata = seu_con$seurat_low.res,
	#weight.reduction = seu_mbv[["pca"]], dims = 1:30)
	weight.reduction = seu_mbv[["PRECAST"]], dims=1:15, k.weight=50)
write.csv(mbv_query, "processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered-conservative_ref-control_query-MBv-conservative_qual-genes-kanchor-50-pc20_red-precast-kweight-50-low-res.csv", row.names=T)
cat("\nSaved label transfer to: processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered-conservative_ref-control_query-MBv-conservative_qual-genes-kanchor-50-pc20_red-precast-kweight-50-low-res.csv\n")

## Reproducibility information
print("\n\nReproducibility information:")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
