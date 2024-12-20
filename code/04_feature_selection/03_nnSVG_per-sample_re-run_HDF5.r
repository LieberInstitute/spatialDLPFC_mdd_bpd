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
set.seed(123)

load(file=paste0("processed-data/04_feature_selection/per-slide_spe/",args[[1]]))
dim(tmp)

cat("\n\nConvert to dgCMatrix...")
format(Sys.time(), tz="EST")
regular_matrix_counts <- as.matrix(assays(tmp)[["counts"]])
sparse_matrix_counts <- as(regular_matrix_counts, "dgCMatrix")
assays(tmp)$counts <- sparse_matrix_counts

cat("\nSubsetting to individual samples...\n")
svg.list = list.files("processed-data/04_feature_selection/per-sample_svgs")
missing.files = setdiff(paste0(tmp$sample_id, "_nnSVG-results.csv"), svg.list)
missing.samples = substr(missing.files, start=0, stop=13)

l1 = unique(tmp$sample_id)
l1 = intersect(l1, missing.samples)
if(length(l1)==0) stop("\nNo missing files\n")
cat("Missing samples:",l1)
names(l1) = l1
spe.list = lapply(l1, function(x) tmp[,tmp$sample_id==x])

cat("\nLooping nnSVG...\n")
for(i in seq_along(spe.list)) {
	cat(names(spe.list)[i])
	cat("\nFilter out zero genes (if present) and generate logcounts on subset of genes...\n")
	#keep_rows = rowSums(counts(spe.list[[i]]))!=0
	#I can't figure out what is causing the error but I know that many of the missing samples have the lowest gene expr
	#so i will see if these complete when I remove any gene that has non-zero expr in <10 spots
	### other runs completed when there were genes with fewer than 10 spots of non-zero expr, but these samples have far more of these genes so maybe that is related??
	#keep_rows = rowSums(counts(spe.list[[i]])!=0)>10
	#tmp_sub = spe.list[[i]][keep_rows,]
	#still wasn't strict enough for some of the slides to run
	tmp_sub = filter_genes(spe.list[[i]], filter_genes_ncounts = 3, filter_genes_pcspots = .1, filter_mito=T)
	cat("\nGenes remaining:",dim(tmp_sub)[1])
	tmp_sub <- computeLibraryFactors(tmp_sub)
	tmp_sub <- logNormCounts(tmp_sub)
	cat(dim(logcounts(tmp_sub)))
	cat("\nnnSVG start...",format(Sys.time(),tz="EST"),"\n")
	set.seed(123)
	results <- nnSVG(tmp_sub, n_threads=12)
	svg = rowData(results)
	cat(names(spe.list)[i],"save output",format(Sys.time(),tz="EST"),"\n\n")
	write.csv(svg, paste0("processed-data/04_feature_selection/per-sample_svgs/", names(spe.list)[i], "_nnSVG-results.csv"), row.names=T)
}



## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
