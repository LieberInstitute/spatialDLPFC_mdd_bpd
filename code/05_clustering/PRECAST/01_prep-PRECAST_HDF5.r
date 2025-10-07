setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(Seurat)
})
set.seed(123)

sampleList = list.files("processed-data/04_feature_selection/per-sample_spe-conservative")
names(sampleList) = substr(sampleList, start=0, stop=13)

#keep.genes = read.csv("processed-data/04_feature_selection/tmp_smaller-feature-list_n1721.csv")
#keep.genes = read.csv("processed-data/04_feature_selection/tmp_larger-feature-list_n3198.csv")
#keep.genes = read.csv("processed-data/04_feature_selection/selected-SVGs_n1051.csv")
geneList = readRDS("processed-data/04_feature_selection/nnSVG-eval_conservative_geneList.rds")
avg.expr = read.csv("processed-data/04_feature_selection/nnSVG-filtered-genes_conservative_avg-logcounts.csv", row.names=1)# %>%
#  tibble::rownames_to_column(var="gene_id")

svg_id = rownames(avg.expr)[avg.expr$gene_name %in% geneList$qual_genes]
length(svg_id)

srt.sets = lapply(sampleList, function(x) {
#for(i in slideList) {
	cat(x,"\n")
	load(paste0("processed-data/04_feature_selection/per-sample_spe-conservative/",x))
	
	#keep.gene.id = rownames(tmp)[rowData(tmp)$gene_name %in% keep.genes$gene_name]
	tmp = tmp[svg_id,]
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
save(srt.sets, file=paste0("processed-data/05_clustering/PRECAST/srt-list_spe-conservative_n",length(svg_id),"_counts.Rdata"))

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
