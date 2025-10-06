args = commandArgs(TRUE)
print(args[[1]])
setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(HDF5Array)
  library(DelayedArray)
  library(scran)
  library(nnSVG)
})
#set.seed(123) #initial run seed == 123
set.seed(456) #re-run seed == 456 to see if that helps with completion of 302_B1 and 333_C1

load(file=paste0("processed-data/04_feature_selection/per-sample_spe-conservative/",args[[1]]))
dim(tmp)

cat("\nFilter out zero genes that are present in <100 spots (per nnSVG documentation recommendation)\nhttps://bioconductor.org/packages/3.21/bioc/vignettes/nnSVG/inst/doc/nnSVG.html#5_Troubleshooting\n")#(if present) and generate logcounts on subset of genes...\n")
keep_rows = rowSums(counts(tmp)>0)>100
table(keep_rows)
tmp_sub = tmp[keep_rows,]
rowData(tmp_sub)$n_spots_nonzero = rowSums(counts(tmp_sub)>0)
#remove any with spots with zero genes now, should be ok to run if only 1 gene
tmp_sub$n_genes_nonzero = colSums(counts(tmp_sub)>0)
tmp_sub = tmp_sub[,tmp_sub$n_genes_nonzero>0]

tmp_sub <- computeLibraryFactors(tmp_sub)
tmp_sub <- logNormCounts(tmp_sub)

cat("\nAll spots have at least 1 non-zero gene\n")
sum(colSums(logcounts(tmp_sub)) > 0)==dim(tmp_sub)[2]
cat("\nMin. number of genes across all spots:",min(colSums(logcounts(tmp_sub)>0)),"\n")
cat("All genes have at least 1 non-zero spot\n")
sum(rowSums(logcounts(tmp_sub)) > 0)==dim(tmp_sub)[1]
cat("\nMin. number of spots across all genes:",min(rowSums(logcounts(tmp_sub)>0)),"\n")

cat("\nStandard nnSVG...",format(Sys.time(),tz="EST"),"\n")
results <- nnSVG(tmp_sub, assay_name="logcounts", n_threads=12)
svg = rowData(results)
cat("\nSave standard nnSVG output...", format(Sys.time(),tz="EST"),"\n")
write.csv(svg, paste0("processed-data/04_feature_selection/per-sample_svgs-conservative/",unique(tmp_sub$sample_id), "_nnSVG-results.csv"), row.names=T)
cat("Saved to:",paste0("processed-data/04_feature_selection/per-sample_svgs-conservative/",unique(tmp_sub$sample_id), "_nnSVG-results.csv"))

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
