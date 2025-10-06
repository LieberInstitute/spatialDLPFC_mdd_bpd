setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(BiocParallel)
	library(nnSVG)
})
set.seed(123)
setAutoBlockSize(1e9)

#cat("\nBiocParallel defaults:\n")
#MulticoreParam()

system.time(spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC-conservative_norm_"))
cat("\nDim spe:",dim(spe),"\n")
cat("\ntree:\n")
showtree(logcounts(spe))
cat("\nseed:\n")
seed(logcounts(spe))

spe = filter_genes(spe, filter_genes_ncounts = 3, filter_genes_pcspots = .1, filter_mito=T) #will give ~6k genes

cat("\nDim spe (filtered genes):",dim(spe),"\n")

#save filtered gene list and decile
avg1 = rowMeans(logcounts(spe))
q.decile = quantile(avg1, prob=seq(0,1,.1))
avg.expr = data.frame("gene_name"=rowData(spe)$gene_name, "avg_expr"=avg1, "decile"=cut(avg1, breaks=c(0,q.decile[2:11]), labels=F))
rownames(avg.expr) = rownames(spe)
write.csv(avg.expr, "processed-data/04_feature_selection/nnSVG-filtered-genes_conservative_avg-logcounts.csv", row.names=T)

#l1 = unique(spe$sample_id)
#names(l1) = l1
#l1 = bplapply(l1, function(x) {
#	print(x); format(Sys.time(), tz="EST")
#	tmp <- spe[,spe$sample_id==x]
#	save(tmp, file=paste0("processed-data/04_feature_selection/per-sample_spe-conservative/",x,".Rdata"))
#})

## Reproducibility information
cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
