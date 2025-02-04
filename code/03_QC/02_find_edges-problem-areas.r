setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(scuttle)
	library(raster)
})
source("code/03_QC/edgeDetection_finalized/raster_edge_functions.r")

system.time(spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_"))
dim(spe)

spe$lg10.umi = log10(spe$sum_umi)
spe$lg10.genes = log10(spe$sum_gene)

cat("\n3MAD outlier to ID poor quality areas:", format(Sys.time()),"\n")
colData(spe)$umi_3MAD.outlier_slide = isOutlier(spe$lg10.umi, subset=spe$in_tissue, batch=spe$slide, type="lower", nmads=3)
colData(spe)$genes_3MAD.outlier_slide = isOutlier(spe$lg10.genes, subset=spe$in_tissue, batch=spe$slide, type="lower", nmads=3)

colData(spe)$umi_3MAD.outlier_sample = isOutlier(spe$lg10.umi, subset=spe$in_tissue, batch=spe$sample_id, type="lower", nmads=3)
colData(spe)$genes_3MAD.outlier_sample = isOutlier(spe$lg10.genes, subset=spe$in_tissue, batch=spe$sample_id, type="lower", nmads=3)

#colData(spe)$umi_3MAD.outlier = ifelse(spe$in_tissue==FALSE, "off tissue", paste(spe$umi_3MAD.outlier_sample, spe$umi_3MAD.outlier_slide))
#colData(spe)$umi_3MAD.outlier = factor(spe$umi_3MAD.outlier,
#                                       levels=c("off tissue","FALSE FALSE","FALSE TRUE","TRUE FALSE","TRUE TRUE"),
#                                       labels=c("off tissue","none","slide","sample","both"))

#colData(spe)$genes_3MAD.outlier = ifelse(spe$in_tissue==FALSE, "off tissue", paste(spe$genes_3MAD.outlier_sample, spe$genes_3MAD.outlier_slide))
#colData(spe)$genes_3MAD.outlier = factor(spe$genes_3MAD.outlier,
#                                         levels=c("off tissue","FALSE FALSE","FALSE TRUE","TRUE FALSE","TRUE TRUE"),
#                                         labels=c("off tissue","none","slide","sample","both"))

colData(spe)$genes_3MAD.outlier_binary = colData(spe)$genes_3MAD.outlier_slide | colData(spe)$genes_3MAD.outlier_sample
colData(spe)$genes_3MAD.outlier_binary = ifelse(colData(spe)$in_tissue==FALSE, FALSE, colData(spe)$genes_3MAD.outlier_binary)

sampleList = unique(spe$sample_id)
names(sampleList) <- sampleList

cat("\nDetect edges:", format(Sys.time()),"\n")
genes_edges = lapply(sampleList, function(x) {
	tmp = colData(spe)[spe$sample_id==x,c("in_tissue","array_row","array_col","genes_3MAD.outlier_binary")]
	clumpEdges(tmp[,-1], rownames(tmp)[tmp$in_tissue==FALSE])
})

cat("\nNumber of samples with edges detected")
table(sapply(genes_edges, length)>0) 
cat("\nSamples with edges detected\n")
names(genes_edges)[sapply(genes_edges, length)>0]
colData(spe)$edge_outlier_genes = FALSE
colData(spe)[unlist(genes_edges),"edge_outlier_genes"] = TRUE

cat("\nFind problem areas:", format(Sys.time()),"\n")
genes_probs = lapply(sampleList, function(x) {
	tmp = colData(spe)[spe$sample_id==x,c("in_tissue","array_row","array_col","genes_3MAD.outlier_binary")]
	problemAreas(tmp[,-1], 
		rownames(tmp)[tmp$in_tissue==FALSE], uniqueIdentifier=x, shifted=F)
})
genes_probs = do.call(rbind, genes_probs)
colData(spe)$problem_areas_genes.id = NA
colData(spe)[genes_probs$spotcode,"problem_areas_genes.id"] = genes_probs$clumpID
colData(spe)$problem_areas_genes.size = 0
colData(spe)[genes_probs$spotcode,"problem_areas_genes.size"] = genes_probs$clumpSize

write.csv(colData(spe), "processed-data/03_QC/colData_edges-problem-areas.csv", row.names=T)
cat("\ncolData saved to: processed-data/03_QC/colData_edges-problem-areas.csv\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()

