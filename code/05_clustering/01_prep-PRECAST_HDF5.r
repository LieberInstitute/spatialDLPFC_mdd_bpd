setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(Seurat)
})
set.seed(123)

slideList = list.files("processed-data/04_feature_selection/per-slide_spe")
srt.sets = list()

for(i in slideList) {
	load(paste0("processed-data/04_feature_selection/per-slide_spe/",i))
	###
	# for purposes of testing limits of PRECAST, just subset to first 2k genes
	tmp = tmp[1:2000,]
	###		
	cat(i,"- Convert counts from HDF5 to dgCMatrix\n")
	regular_matrix_counts <- as.matrix(assays(tmp)[["counts"]])
	sparse_matrix_counts <- as(regular_matrix_counts, "dgCMatrix")
	assays(tmp)$counts <- sparse_matrix_counts
	
	l2 = unique(tmp$sample_id)
	names(l2) = lapply(l2, function(x) unique(colData(tmp)[tmp$sample_id==x,"brain"]))
	l2 = lapply(l2, function(x) tmp[,colData(tmp)$sample_id==x])

	cat(i,"- per-sample SPE to per-sample SeuratObject\n")
	srt.sets = c(srt.sets, lapply(l2, function(x) {
		rownames(colData(x)) <- paste(x$sample_id, rownames(colData(x)), sep="_")
		colnames(counts(x)) <- rownames(colData(x))
		colData(x)$col <- x$array_col
		colData(x)$row <- x$array_row
		count <- counts(x)
		a1 <- CreateAssayObject(count, assay = "RNA", min.features = 0, min.cells = 0)
		CreateSeuratObject(a1, meta.data = as.data.frame(colData(x)))
	}))
}
cat("\n\nFinal srt.sets structure:\n")
str(srt.sets, 3)
save(srt.sets, file="processed-data/05_clustering/srt-list_spe-HDF5_counts.Rdata")

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
