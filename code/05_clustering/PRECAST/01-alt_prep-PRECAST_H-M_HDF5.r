setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(Seurat)
	library(dplyr)
})
set.seed(123)

sampleList = list.files("processed-data/04_feature_selection/per-sample_spe")
names(sampleList) = substr(sampleList, start=0, stop=13)

layer.markers = read.csv("processed-data/04_feature_selection/EXT_TableS9_sig_genes_FDR5perc_enrichment.csv") %>%
	filter(stat>0, spatial_domain_resolution=="Sp09") %>%
	mutate(domain_simple=factor(test, 
		levels=paste0("Sp09D0",c(1,2,3,5,8,4,7,6,9)), 
		labels=c("L1 (1)","L1 (2)","L2","L3","L4","L5","L6","WM","WM")))
domains = c("L1 (1)","L1 (2)","L2","L3","L4","L5","L6","WM")
names(domains) = domains
top500 = lapply(domains, function(x) filter(layer.markers, domain_simple==x) %>% slice_min(n=500, fdr) %>% pull(gene))
cat("\nNumber of unique genes from the top 500 lowest FDR per Huuki-Myers Sp09 domain:\n")
length(unique(unlist(top500)))

avg.expr = read.csv("processed-data/04_feature_selection/nnSVG-filtered-genes_avg-logcounts.csv", row.names=1)# %>%
#  tibble::rownames_to_column(var="gene_id")
top500 = lapply(top500, intersect, y=avg.expr$gene_name)
cat("\nNumber of top 500 layer markers (non-unique) also in 6111 pre-filtered genes:\n")
sapply(top500, length)
top500_names = unique(unlist(top500))
top500_ids = rownames(avg.expr)[avg.expr$gene_name %in% top500_names]
cat("\nNumber of unique top 500 layer markers also in 6111 pre-filtered genes:\n")
length(top500_ids)

srt.sets = lapply(sampleList, function(x) {
#for(i in slideList) {
	cat(x,"\n")
	load(paste0("processed-data/04_feature_selection/per-sample_spe/",x))
	
	#keep.gene.id = rownames(tmp)[rowData(tmp)$gene_name %in% keep.genes$gene_name]
	tmp = tmp[top500_ids,]
	dim(tmp)

	regular_matrix_counts <- as.matrix(assays(tmp)[["counts"]])
	sparse_matrix_counts <- as(regular_matrix_counts, "dgCMatrix")
	assays(tmp)$counts <- sparse_matrix_counts
	
	#l2 = unique(tmp$sample_id)
	#names(l2) = lapply(l2, function(x) unique(colData(tmp)[tmp$sample_id==x,"brain"]))
	#l2 = lapply(l2, function(x) tmp[,colData(tmp)$sample_id==x])

#	cat(i,"- per-sample SPE to per-sample SeuratObject\n")
#	srt.sets = c(srt.sets, lapply(l2, function(x) {
		rownames(colData(tmp)) <- paste(tmp$sample_id, rownames(colData(tmp)), sep="_")
		colnames(counts(tmp)) <- rownames(colData(tmp))
		colData(tmp)$col <- tmp$array_col
		colData(tmp)$row <- tmp$array_row
		count <- counts(tmp)
		a1 <- CreateAssayObject(count, assay = "RNA", min.features = 0, min.cells = 0)
		CreateSeuratObject(a1, meta.data = as.data.frame(colData(tmp)))
#	}))
})
cat("\n\nFinal srt.sets structure:\n")
str(srt.sets, 3)
save(srt.sets, file="processed-data/05_clustering/PRECAST/srt-list_spe_H-M-markers_counts.Rdata")

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
