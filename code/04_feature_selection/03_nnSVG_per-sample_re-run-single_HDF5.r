setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(HDF5Array)
  library(DelayedArray)
  library(scran)
  library(nnSVG)
})
set.seed(123)

load(file="processed-data/04_feature_selection/per-slide_spe/V13B23-331.Rdata")
dim(tmp)

cat("\n\nConvert to dgCMatrix...")
format(Sys.time(), tz="EST")
regular_matrix_counts <- as.matrix(assays(tmp)[["counts"]])
sparse_matrix_counts <- as(regular_matrix_counts, "dgCMatrix")
assays(tmp)$counts <- sparse_matrix_counts

cat("\nSubsetting to sample V13B23-331_D1...\n")
tmp_sub = tmp[,tmp$sample_id=="V13B23-331_D1"]

cat("\nFilter out zero genes (if present) and generate logcounts on subset of genes...\n")
tmp_sub = filter_genes(tmp_sub, filter_genes_ncounts = 3, filter_genes_pcspots = .1, filter_mito=T)
tmp_sub <- computeLibraryFactors(tmp_sub)
tmp_sub = tmp_sub[,sizeFactors(tmp_sub)>0]
tmp_sub <- logNormCounts(tmp_sub)
cat(dim(logcounts(tmp_sub)))

cat("\nnnSVG start...",format(Sys.time(),tz="EST"),"\n")
set.seed(123)
results <- nnSVG(tmp_sub, n_threads=12)
svg = rowData(results)
cat("V13B23-331_D1 save output",format(Sys.time(),tz="EST"),"\n\n")
write.csv(svg, "processed-data/04_feature_selection/per-sample_svgs/V13B23-331_D1_nnSVG-results.csv", row.names=T)



## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
