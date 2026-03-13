suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(HDF5Array)
  library(ggplot2)
  library(dplyr)
  library(loomR)
})
set.seed(123)
setAutoBlockSize(1e9)

#load spe
cat("\nFull spe...\n")
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")
dim(spe) # 28965 535248


#load clusters
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
stopifnot(identical(rownames(colData(spe)), rownames(cdata)))
cdata$smoothed_k9_1663 = factor(cdata$smoothed_k9_1663_f, levels=c("L1","L2","L3/4","L5","L6","WM","low UMI","Vasc","GABA"),
                                labels=c("L1","L2","L3.4","L5","L6","WM","low.UMI","Vasc","GABA"))
cdata$key = colData(spe)$key
cdata$nUMI = colData(spe)$sum_umi
cdata$nGene = colData(spe)$sum_gene

# add seurat label
res = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc30_red-precast-kweight-50-low-res.csv", row.names=1)
res = res[rownames(cdata),]
stopifnot(identical(rownames(res), rownames(cdata)))
cdata$seurat_label = factor(res$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"),
                            labels=c("Micro.Vasc","Astro","L2.3","L2.3","L4","Inhb","L5","L6","Oligo"))

#for rowdata want to import pseudobulk edgeR filterByExpr results
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-no-lowUMI_sample-seurat-pc30_norm-filt.Rdata")
cat("\nPseudobulk spe for gene filtering...\n")
dim(spe_pseudo) # 21080   690
rowData(spe_pseudo)$high_expr_group_sample_id2 <- edgeR::filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster2 <- edgeR::filterByExpr(spe_pseudo, group = spe_pseudo$seurat_label)

rdata = as.data.frame(rowData(spe_pseudo)[,c("gene_id","gene_name","gene_type")])
rdata$DE_input = rowData(spe_pseudo)$high_expr_group_sample_id2==T & rowData(spe_pseudo)$high_expr_group_cluster2==T
table(rdata$DE_input)

cat("\nFiltered gene set to DE input...\n")
rdata = rdata[rdata$DE_input,]
dim(rdata)

spe = spe[rownames(rdata),]

# correct replicated gene names
rowData(spe)[c("ENSG00000187522"),"gene_name"] = "MSTANTD7"
rdata2 = rowData(spe)[,c("gene_id","gene_name","gene_type")]
rdata2$DE_input = rdata[rownames(rdata2),"DE_input"]
rownames(rdata2) = rdata2$gene_name
cat("\nNumber of genes after revising to unique gene names:\n")
nrow(rdata2)


# remove low UMI cluster
cat("\nRemoving low UMI cluster...\n")
cdata = cdata[cdata$smoothed_k9_1663!="low.UMI",]
cdata$smoothed_k9_1663 = droplevels(cdata$smoothed_k9_1663)

spe = spe[,rownames(cdata)]
dim(spe)
# 513200 spots

if(file.exists("processed-data/09_SCENIC/cdata_revised-qc-metrics_for-AUCell-13162.csv")) {
	out1 <- read.csv("processed-data/09_SCENIC/cdata_revised-qc-metrics_for-AUCell-13162.csv", row.names=1)
	out1$smoothed_k9_1663 = factor(out1$smoothed_k9_1663, levels=levels(cdata$smoothed_k9_1663))
	out1$seurat_label = factor(out1$seurat_label, levels=levels(cdata$seurat_label))
} else {
	out1 = scuttle::perCellQCMetrics(spe, assay.type="counts")
	stopifnot(identical(rownames(cdata), rownames(out1)))
	out1 = cbind(out1, cdata[,c("sex","condition","smoothed_k9_1663","seurat_label")])
	write.csv(out1, "processed-data/09_SCENIC/cdata_revised-qc-metrics_for-AUCell-13162.csv", row.names=T)
}

table(out1$smoothed_k9_1663)
cat("\n")
sapply(levels(out1$smoothed_k9_1663), function(x) quantile(filter(as.data.frame(out1), smoothed_k9_1663==x)$detected, probs=c(.01,.05,.1,.5,1)))
cat("\n\n")

table(out1$seurat_label)
cat("\n")
sapply(levels(out1$seurat_label), function(x) quantile(filter(as.data.frame(out1), seurat_label==x)$detected, probs=c(.01,.05,.1,.5,1)))
cat("\n\n")


colnames(out1)[grep("sum", colnames(out1))] = "nUMI"
colnames(out1)[grep("detected", colnames(out1))] = "nGene"

#further split astro and L23 
## splitting astro for better sampled representation of # detected genes 
## splitting L23 because they are equal in # detected genes and cutting down size will make for more efficient array job memory requests
#out1$custom_cluster = as.character(out1$seurat_label)
#out1[out1$custom_cluster=="Astro" & out1$smoothed_k9_1663 %in% c("L1","Vasc","WM"), "custom_cluster"] = "Astro.Glia"
#out1[out1$custom_cluster=="Astro" & !out1$smoothed_k9_1663 %in% c("L1","Vasc","WM"), "custom_cluster"] = "Astro.Nrn"
#out1[out1$custom_cluster=="L2.3" & out1$smoothed_k9_1663=="L2", "custom_cluster"] = "L2"
#out1[out1$custom_cluster=="L2.3" & out1$smoothed_k9_1663!="L2", "custom_cluster"] = "L3"
#
#out1$custom_cluster = factor(out1$custom_cluster, levels=c("Micro.Vasc","Astro.Glia","Astro.Nrn","L2","L3","L4","Inhb","L5","L6","Oligo"))
#
#table(out1$custom_cluster)
#cat("\n")
#sapply(levels(out1$custom_cluster), function(x) quantile(filter(as.data.frame(out1), custom_cluster==x)$nGene, probs=c(.01,.05,.1,.5,1)))
#cat("\n\n")
#
#write.csv(out1, "processed-data/09_SCENIC/cdata_revised-qc-metrics_custom-cluster_for-AUCell.csv", row.names=T)

#for(i in c("Astro.Glia","Astro.Nrn","L2","L3")) {
for(i in levels(out1$seurat_label)) {

	cat(paste0("\n\n",i,"...\n"))

#        (fn = paste0("processed-data/09_SCENIC/expr_loom/spe-n119_custom-cluster-",i,"_", nrow(rdata2), "-genes_no-lowUMI_logcounts.loom"))
	(fn = paste0("processed-data/09_SCENIC/expr_loom/spe-n119_seurat-label-",i,"_", nrow(rdata2), "-genes_no-lowUMI_logcounts.loom"))

#        tmp = out1[out1$custom_cluster==i,]
	tmp = out1[out1$seurat_label==i,]
	nrow(tmp)

	mtx = logcounts(spe[,rownames(tmp)])
	stopifnot(identical(rownames(mtx), rdata2$gene_id))
	rownames(mtx) = rdata2$gene_name
	dim(mtx)

	lfile = create(fn, mtx,
		cell.attrs=as.list(tmp[colnames(mtx),c("condition","sex","smoothed_k9_1663","seurat_label","nUMI","nGene")]),
		calc.numi = F, chunk.dims=c(100,500))
	stopifnot(identical(rownames(rdata2), lfile[["row_attrs/Gene"]][]))
	lfile$add.row.attribute(list("gene_id" = rdata2$gene_id,
                             "gene_name" = rownames(rdata2),
                             "gene_type" = rdata2$gene_type,
                             "DE_input" = rdata2$DE_input), 
                        overwrite = TRUE)
	lfile$close_all()
	cat(paste0("\n",i," saved!\n"))
}

## Reproducibility information
cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
