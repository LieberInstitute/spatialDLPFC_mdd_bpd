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

#all candidate SVGs
id.1663 = rownames(spe)[rowData(spe)$gene_name %in% geneList$qual_genes]
cat("\n\nQualifying genes (including batch effect and low spcov):\n")
length(unique(id.1663))

format(Sys.time(), tz="EST")
spe <- runPCA(spe, subset_row=id.1663, ncomponents=50, name="PCA_1663", BSPARAM=RandomParam(), BPPARAM=MulticoreParam(workers=12))
reducedDimNames(spe)

#remove top decile low spcov
id.1104 = rownames(spe)[rowData(spe)$gene_name %in% setdiff(geneList$qual_genes, geneList$top.decile_low.spcov)]
cat("\n\nQualifying genes (top decile low spcov):\n")
length(unique(id.1104))

format(Sys.time(), tz="EST")
set.seed(123)
spe <- runPCA(spe, subset_row=id.1104, ncomponents=50, name="PCA_1104", BSPARAM=RandomParam(), BPPARAM=MulticoreParam(workers=12))
reducedDimNames(spe)

#remove batch effect genes
id.1079	= rownames(spe)[rowData(spe)$gene_name %in% geneList$final_svgs]
cat("\n\nQualifying genes (without batch effect	genes, low spcov high expr. genes removed):\n")
length(unique(id.1079))

format(Sys.time(), tz="EST")
set.seed(123)
spe <- runPCA(spe, subset_row=id.1079, ncomponents=50, name="PCA_1079", BSPARAM=RandomParam(), BPPARAM=MulticoreParam(workers=12))
reducedDimNames(spe)

#dummy_spe = spe
#assay(dummy_spe, "counts") <- NULL
#assay(dummy_spe, "logcounts") <- NULL
#dummy_spe
#cat("\n\nSave SpatialExperiment with empty assays (to access reducedDims)\n")
#format(Sys.time(), tz="EST")
#save(dummy_spe, file="processed-data/04_feature_selection/spe_PCA_empty-assays.Rdata")
#cat("\nSaved to: processed-data/04_feature_selection/spe_PCA_empty-assays.Rdata")
quickResaveHDF5SummarizedExperiment(spe)
cat("\nSaved with quickResave to spe HDF5 dir\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
