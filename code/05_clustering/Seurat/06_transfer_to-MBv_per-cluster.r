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

#find anchors per cluster
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
stopifnot(identical(rownames(seu_mbv@meta.data), rownames(cdata)))
seu_mbv@meta.data$precast_k9_1663_f = factor(cdata$precast_k9_1663_f, levels=c("Vasc","L1","L2","L3","GABA","L5","L6","WM","low UMI"))
Idents(seu_mbv) <- "precast_k9_1663_f"
seu_list = SplitObject(seu_mbv, split.by="ident")
for(i in levels(seu_mbv@meta.data$precast_k9_1663_f)) {
	cat("\n\nProcessing", i, "...\n")
	tmp = seu_list[[i]]
	tmp = SCTransform(tmp, verbose=T)
	cat("\nStart integration...\n")
	format(Sys.time())
	anchors <- FindTransferAnchors(reference = seu_con, query = tmp, normalization.method = "SCT",
		k.anchor=20, k.score=50, max.features=500, npcs=50, dims=1:50, mapping.score.k=T)
	save(anchors, file=paste0("processed-data/05_clustering/Seurat/anchors_SZBDMulti-seq_MBv-filtered_ref-control_query-MBv-",gsub(" ","-",i),".Rdata"))
	cat("\nSaved to:",paste0("processed-data/05_clustering/Seurat/anchors_SZBDMulti-seq_MBv-filtered_ref-control_query-MBv-",gsub(" ","-",i),".Rdata"),"\n")
}

##load anchors
#cat("\nLoad saved anchors...\n")
#load("processed-data/05_clustering/Seurat/anchors_SZBDMulti-seq_MBv-filtered_ref-control_query-MBv-kanchor-50.Rdata")

#fix seu_con cluster resolution
#cat("\nRe-do seu_con cluster resolution and annotate...\n")
#format(Sys.time())
#set.seed(123)
#seu_con <- FindClusters(seu_con, resolution=.3) 

#azmap = c("L3"="0","L2"="1","L5"="2","Oligo"="3","L4/L5"="4","L6"="5",
#          "Inhb CGE VIP"="6","Inhb MGE PV"="7", "L6b"="8","Inhb MGE SST"="9","Astro"="10",
#          "L3/L4"="11","Inhb CGE LAMP5"="12","OPC"="13","L5/6 NP"="14","Micro/Vasc"="15",
#          "L6 IT Car3"="16", "Chandelier"="17", "L5 ET"="18")
##azmap = c("L2/3 A"="0","L5 IT"="1","Inhb A"="2","Oligo"="3","L2/3 B"="4","L4 A"="5",
##          "L6 IT"="6","Inhb B"="7", "L6 CT/ L6b"="8","Astro/Vasc"="9","L4 B"="10",
##          "Inhb LAMP5"="11","OPC"="12","L5/6 NP"="13","L6 IT Car3"="14","Micro/Immune"="15",
##          "Chandelier"="16", "L5 ET"="17")
#seu_con$seurat_annotated = factor(as.character(seu_con$seurat_clusters), levels=azmap,
#                                  labels=names(azmap))

##seurat_high.res = factor(as.character(seu_con$seurat_annotated),
##                                 levels=c("Astro/Vasc","Micro/Immune",
##                                          "Inhb A","Chandelier","Inhb B","Inhb LAMP5",
##                                          "L2/3 B","L2/3 A",
##                                          "L4 B","L4 A",
##                                          "L5 IT","L5 ET","L5/6 NP",
##                                          "L6 IT","L6 IT Car3",
##                                          "L6 CT/ L6b",
##                                          "Oligo","OPC"),
##                                 labels=c("Astro/Vasc","Micro/Immune",
##                                          "Inhb MGE","Inhb MGE","Inhb CGE","Inhb LAMP5",
##                                          "L2","L3",
##                                          "L3/L4","L4/L5",
##                                          "L5", "L5", "L5",
##                                          "L6", "L6", 
##                                          "L6b",
##                                          "Oligo", "Oligo"))

#seurat_low.res <- factor(as.character(seu_con$seurat_annotated),
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

##seu_con$seurat_low.res = factor(as.character(seu_con$seurat_high.res),
##                                levels=c("Astro/Vasc","Micro/Immune",
##                                         "Inhb MGE","Inhb CGE","Inhb LAMP5",
##                                         "L2","L3",
##                                         "L3/L4","L4/L5",
##                                         "L5", "L6", "L6b",
##                                         "Oligo"),
##                                labels=c("Glia (non-olig)","Glia (non-olig)",
##                                         "Inhb","Inhb","Inhb",
##                                         "L2","L3",
##                                         "L4","L4",
##                                         "L5","L6","L6b",
##                                         "Oligo"))

#rm(seu_con) #to save space for global vars (not sure if helps but can't hurt)

#cat("\nLoad in PRECAST results for query weights...\n")
#format(Sys.time())
#load("processed-data/05_clustering/PRECAST/srt_precast_k-9_n1663.Rdata")
##precast colnames are different and reordered
#emb_mtx = seuInt[["PRECAST"]]@cell.embeddings
#seurat_key = paste(seu_mbv$sample_id, colnames(seu_mbv), sep="_")
#emb_mtx = emb_mtx[seurat_key,]
#rownames(emb_mtx) = colnames(seu_mbv)

#cat("\nTransfer labels...\n")
#format(Sys.time())

#mbv_query <- TransferData(anchorset = anchors, refdata = seurat_low.res, #prediction.assay = TRUE,
#	#weight.reduction = seu_mbv[["pca"]], dims = 1:30)
#	weight.reduction = emb_mtx, dims=1:15)
#write.csv(mbv_query, "processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_kanchor-50-PRECAST-weights_seurat-low-res.csv", row.names=T)
#cat("\nSaved label transfer to: processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_kanchor-50-PRECAST-weights_seurat-low-res.csv\n")

## Reproducibility information
print("\n\nReproducibility information:")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
