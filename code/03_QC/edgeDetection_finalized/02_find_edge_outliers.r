setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(raster)
})
source("code/03_QC/edgeDetection_finalized/raster_edge_functions.r")

cdata = read.csv("processed-data/03_QC/edgeDetection_finalized/colData_3MAD.csv", row.names=1)
cdata$umi_3MAD.outlier_binary = cdata$umi_3MAD.outlier_slide | cdata$umi_3MAD.outlier_sample
cdata$umi_3MAD.outlier_binary = ifelse(cdata$in_tissue==FALSE, FALSE, cdata$umi_3MAD.outlier_binary)
cdata$genes_3MAD.outlier_binary = cdata$genes_3MAD.outlier_slide | cdata$genes_3MAD.outlier_sample
cdata$genes_3MAD.outlier_binary = ifelse(cdata$in_tissue==FALSE, FALSE, cdata$genes_3MAD.outlier_binary)
table(cdata[,c("umi_3MAD.outlier_binary","genes_3MAD.outlier_binary")])

sampleList = unique(cdata$sample_id)
names(sampleList) <- sampleList

umi_edges = lapply(sampleList, function(x) {
  clumpEdges(cdata[cdata$sample_id==x,c("array_row","array_col","umi_3MAD.outlier_binary")])
})
length(umi_edges) #120
table(sapply(umi_edges, length)>0) #32 samples with detected edges for umi outliers
cdata$edge_outlier_umi = FALSE
cdata[unlist(umi_edges),"edge_outlier_umi"] = TRUE

genes_edges = lapply(sampleList, function(x) {
  clumpEdges(cdata[cdata$sample_id==x,c("array_row","array_col","genes_3MAD.outlier_binary")])
})
length(genes_edges) #120
table(sapply(genes_edges, length)>0) #33 samples with detected edges for genes outliers
cdata$edge_outlier_genes = FALSE
cdata[unlist(genes_edges),"edge_outlier_genes"] = TRUE

genes_edges = lapply(sampleList, function(x) {
  clumpEdges(cdata[cdata$sample_id==x,c("array_row","array_col","genes_3MAD.outlier_binary")], shifted=TRUE)
})
length(genes_edges) #120
table(sapply(genes_edges, length)>0) #33 samples with detected edges for genes outliers
cdata$edge_outlier_genes.shifted = FALSE
cdata[unlist(genes_edges),"edge_outlier_genes.shifted"] = TRUE

write.csv(cdata, "processed-data/03_QC/edgeDetection_finalized/colData_found-edges.csv", row.names=T)
cat("\nNew colData csv with edge outliers saved to: processed-data/03_QC/edgeDetection_finalized/colData_found-edges.csv\n")

cat("\nComparison of UMI-based and # genes-based edge detection\n(# of spots identified)\n")
table(cdata[,c("edge_outlier_umi","edge_outlier_genes")])

cat("\nSample IDs where UMI edges > # genes edges\n(# of spots identified)\n")
cd1 = cdata[cdata$edge_outlier_umi==TRUE & cdata$edge_outlier_genes==FALSE,]
table(cd1$sample_id)
cat("\nSample IDs where UMI edge < # genes edges\n(# of spots identified)\n")
cd2 = cdata[cdata$edge_outlier_umi==FALSE & cdata$edge_outlier_genes==TRUE,]
table(cd2$sample_id)



cat("\nComparison of unshifted (default) and shifted # genes-based edge detection\n(# of spots identified)\n")
table(cdata[,c("edge_outlier_genes","edge_outlier_genes.shifted")])

cat("\nSample IDs where shifted > unshifted # genes edges\n(# of spots identified)\n")
cd1 = cdata[cdata$edge_outlier_genes.shifted==TRUE & cdata$edge_outlier_genes==FALSE,]
if(nrow(cd1)>0) table(cd1$sample_id)
cat("\nSample IDs where shifted < unshifted # genes edges\n(# of spots identified)\n")
cd2 = cdata[cdata$edge_outlier_genes.shifted==FALSE & cdata$edge_outlier_genes==TRUE,]
if(nrow(cd2)>0) table(cd2$sample_id)


## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()

