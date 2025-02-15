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

load(file=paste0("processed-data/04_feature_selection/per-sample_spe/",args[[1]]))
dim(tmp)

cat("\n\nConvert to dgCMatrix...")
format(Sys.time(), tz="EST")
regular_matrix_counts <- as.matrix(assays(tmp)[["counts"]])
sparse_matrix_counts <- as(regular_matrix_counts, "dgCMatrix")
assays(tmp)$counts <- sparse_matrix_counts

#cat("\nSubsetting to individual samples...\n")
#format(Sys.time(), tz="EST")
#l1 = unique(tmp$sample_id)
#names(l1) = l1
#spe.list = lapply(l1, function(x) tmp[,tmp$sample_id==x])

#cat("\nLooping nnSVG...\n")
#for(i in seq_along(spe.list)) {
	#cat(names(spe.list)[i])
	cat("\nFilter out zero genes that are present in <100 spots (per nnSVG documentation recommendation)\nhttps://bioconductor.org/packages/3.21/bioc/vignettes/nnSVG/inst/doc/nnSVG.html#5_Troubleshooting\n")#(if present) and generate logcounts on subset of genes...\n")
	keep_rows = rowSums(counts(tmp)>0)>100
	table(keep_rows)
	tmp_sub = tmp[keep_rows,]
	tmp_sub <- computeLibraryFactors(tmp_sub)
	tmp_sub <- logNormCounts(tmp_sub)
	rowData(tmp_sub)$n_spots_nonzero = rowSums(counts(tmp_sub)>0)

	cat("\nAll spots have at least 1 non-zero gene\n")
	sum(colSums(logcounts(tmp_sub)) > 0)==dim(tmp_sub)[2]
	cat("\nMin. number of genes across all spots:",min(colSums(logcounts(tmp_sub)>0)),"\n")
	cat("All genes have at least 1 non-zero spot\n")
	sum(rowSums(logcounts(tmp_sub)) > 0)==dim(tmp_sub)[1]
	cat("\nMin. number of spots across all genes:",min(rowSums(logcounts(tmp_sub)>0)),"\n")
	
#cat("\ngenerate_weights...",format(Sys.time(),tz="EST"),"\n")
#weights <- generate_weights(input = tmp_sub, stabilize = TRUE, n_threads=12)
#cat("\nSave weights...", format(Sys.time(), tz="EST"),"\n")
#saveRDS(weights, paste0("processed-data/04_feature_selection/per-sample_weights/",unique(tmp_sub$sample_id),"_spoon-weights.rda"))
#format(Sys.time(), tz="EST")

#if(class(logcounts(tmp_sub))[1]=="DelayedMatrix") weights <- DelayedArray(weights)
#source("code/04_feature_selection/spoon-tutorial_test/weighted-nnSVG_re-write_iterative-cov-matrix.r")
#weighted_logcounts <- t(weights)*assays(tmp_sub)[['logcounts']]
#assay(tmp_sub, "weighted_logcounts") <- weighted_logcounts

	cat("\nStandard nnSVG...",format(Sys.time(),tz="EST"),"\n")
	set.seed(123)
	results <- nnSVG(tmp_sub, assay_name="logcounts", n_threads=12)
	svg = rowData(results)
	cat("\nSave standard nnSVG output...", format(Sys.time(),tz="EST"),"\n")
	write.csv(svg, paste0("processed-data/04_feature_selection/per-sample_svgs/",unique(tmp_sub$sample_id), "_nnSVG-results.csv"), row.names=T)
	cat("Saved to:",paste0("processed-data/04_feature_selection/per-sample_svgs/",unique(tmp_sub$sample_id), "_nnSVG-results.csv"))
#}



## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
