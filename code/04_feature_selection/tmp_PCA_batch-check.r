setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(scater)
	library(BiocSingular)
	library(BiocParallel)
})
set.seed(123)

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n120_postQC_norm_")

geneList = readRDS("processed-data/04_feature_selection/nnSVG-eval_geneList.rds")

#top 1600 genes with batch genes included
top.1638 = union(geneList$qual_genes, geneList$qual_batch_effect)
id.1638 = rownames(spe)[rowData(spe)$gene_name %in% top.1638]
cat("\n\nQualifying genes (including batch effect genes):\n")
length(unique(id.1638))

format(Sys.time(), tz="EST")
spe <- runPCA(spe, subset_row=id.1638, ncomponents=50, name="PCA_1638", BSPARAM=RandomParam(), BPPARAM=MulticoreParam(workers=12))
reducedDimNames(spe)

#top 1600 genes without batch effect
top.1604 = geneList$qual_genes
id.1604 = rownames(spe)[rowData(spe)$gene_name %in% top.1604]
cat("\n\nQualifying genes (without batch effect genes):\n")
length(unique(id.1604))

format(Sys.time(), tz="EST")
set.seed(123)
spe <- runPCA(spe, subset_row=id.1604, ncomponents=50, name="PCA_1604", BSPARAM=RandomParam(), BPPARAM=MulticoreParam(workers=12))
reducedDimNames(spe)

#remove highly most expressed genes
top.1054 = setdiff(geneList$qual_genes, geneList$top.decile_low.spcov)
id.1054	= rownames(spe)[rowData(spe)$gene_name %in% top.1054]
cat("\n\nQualifying genes (without batch effect	genes, low spcov high expr. genes removed):\n")
length(unique(id.1054))

format(Sys.time(), tz="EST")
set.seed(123)
spe <- runPCA(spe, subset_row=id.1054, ncomponents=50, name="PCA_1054", BSPARAM=RandomParam(), BPPARAM=MulticoreParam(workers=12))
reducedDimNames(spe)

dummy_spe = spe
assay(dummy_spe, "counts") <- NULL
assay(dummy_spe, "logcounts") <- NULL
dummy_spe
cat("\n\nSave SpatialExperiment with empty assays (to access reducedDims)\n")
format(Sys.time(), tz="EST")
save(dummy_spe, file="processed-data/04_feature_selection/spe_PCA_empty-assays.Rdata")
cat("\nSaved to: processed-data/04_feature_selection/spe_PCA_empty-assays.Rdata")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
