setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(Seurat)
})
set.seed(123)
setAutoBlockSize(1e9)

#spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC-conservative_norm_")
#cat("\nDim spe:",dim(spe),"\n")
#cat("\ntree:\n")
#showtree(logcounts(spe))
#cat("\nseed:\n")
#seed(logcounts(spe))

##realize counts
#cat("\nRealize counts matrix...\n")
#format(Sys.time())
#regular_matrix_counts <- as.matrix(assays(spe)[["counts"]])
#sparse_matrix_counts <- as(regular_matrix_counts, "dgCMatrix")

##subset to shared genes and use gene name
#cat("\nSubset counts to pseudobulk filtered genes...\n")
#load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_processed-SCT.Rdata")

#rdata = read.csv("processed-data/06_pseudobulk/PRECAST_smoothed/pseudobulk-sample-smoothed-n1626-k7_conservative_filtered-genes_avg-logcounts.csv", row.names=1)
#head(rdata)

#fdata = seu_con[["RNA"]]@meta.data
#rownames(fdata) <- rownames(seu_con[["RNA"]])
#head(fdata)

#fdata = fdata[fdata$featureid %in% intersect(fdata$featureid, rownames(rdata)),]

#rdata = rdata[fdata$featureid,]
#dim(rdata)
#stopifnot(identical(rownames(rdata), fdata$featureid))

#rdata$featurekey = rownames(fdata)

##subset counts matrix
#sparse_matrix_counts <- sparse_matrix_counts[rownames(rdata),]
#rownames(sparse_matrix_counts) <- rdata$featurekey
#dim(sparse_matrix_counts)

##finish rdata for seurat mbv feature metadata
#rdata$gene_id = rownames(rdata)
#rownames(rdata) = rdata$featurekey

##create seurat
#cat("\nCreate seurat object...\n")
#format(Sys.time())
#seu_mbv <- CreateSeuratObject(counts = sparse_matrix_counts, meta.data = as.data.frame(colData(spe)))
#seu_mbv[["RNA"]] <- AddMetaData(object = seu_mbv[["RNA"]], metadata = rdata)
#seu_mbv
#save(seu_mbv, file="processed-data/05_clustering/Seurat/seurat_MBv_conservative.Rdata")
#cat("\nFile save checkpoint: processed-data/05_clustering/Seurat/seurat_MBv_conservative.Rdata\n")

#load saved raw seurat object
cat("\nLoading saved, unprocessed seurat object: processed-data/05_clustering/Seurat/seurat_MBv_conservative.Rdata\n")
load("processed-data/05_clustering/Seurat/seurat_MBv_conservative.Rdata")
dim(seu_mbv)

# need to subset genes
load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_conservative_pseudo_sample-smoothed-n1626-k7_norm-filt.Rdata")
pass.genes = rowSums(as.data.frame(rowData(spe_pseudo)[,c("high_expr_group_sample_id","high_expr_group_cluster","high_expr_SZBDMultiseq")]))>2
seu_mbv = seu_mbv[rowData(spe_pseudo)$gene_name[pass.genes],]
cat("\nSubset to genes that passed edgeR filters for PRECAST (smoothed) pseudobulk and are present in SZBDMulti-seq:\n")
dim(seu_mbv)

# need to exclude any spots with no genes now
seu_mbv$nFeature_RNA = colSums(seu_mbv@assays$RNA$counts>0)
seu_mbv = seu_mbv[,seu_mbv$nFeature_RNA>0]
cat("\nRemove any spots that have no counts with filtered genes:\n")
dim(seu_mbv)

#pre-process SRT
seu_mbv <- SCTransform(seu_mbv, ncells=10000, verbose=T) #default ncells is 5k
### ncells parameter is the number of cells on which the noise is learned
### considering that I have 500k, maybe it will help to improve the ncells
seu_mbv <- RunPCA(seu_mbv)
seu_mbv
save(seu_mbv, file="processed-data/05_clustering/Seurat/seurat_MBv_conservative_processed-SCT.Rdata")
cat("\nSaved to: processed-data/05_clustering/Seurat/seurat_MBv_conservative_processed-SCT.Rdata\n")

## Reproducibility information
print("Reproducibility information:")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
