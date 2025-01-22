setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(STexampleData)
	library(scran)
	library(nnSVG)
	library(spoon)
})
set.seed(123)

spe <- Visium_mouseCoronal()
spe <- spe[, colData(spe)$in_tissue == 1]
spe <- filter_genes(spe)
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)

rowData(spe)$avg.norm.expr = rowMeans(logcounts(spe))
rowData(spe)$decile = cut(rowData(spe)$avg.norm.expr, breaks=c(0,quantile(rowData(spe)$avg.norm.expr, probs=seq(0,1,by=.1))[2:11]), labels=F)

nList = c(50,100,500,1000,5000)
names(nList) <- paste0("n",nList)

cat("\nsampleList",format(Sys.time(), tz="EST"),"\n")
sampleList <- list()

for (i in seq_along(nList)) {
  if(i==1) {
    #random sample genes
    sampleList[[i]] <- sample(rownames(spe), nList[[i]])
    #name sampleList level
    names(sampleList)[i] = paste0("n", nList[[i]])
  } else {
    #number of new random genes to sample
    n_random = nList[[i]]-nList[[(i-1)]]
    print(n_random)
    #remaining spe genes eligible to sample
    genes_eligible = setdiff(rownames(spe), sampleList[[(i-1)]])
    #addition random sample of remaining genes plus existing genes 
    sampleList[[i]] <- c(sampleList[[(i-1)]], sample(genes_eligible, n_random))
    #name sampleList level
    names(sampleList)[i] = paste0("n", nList[[i]])
  } 
}
sapply(sampleList, length)
table(table(unlist(sampleList)))

cat("\nspeList",format(Sys.time(), tz="EST"),"\n")
speList <- lapply(sampleList, function(x) {
  tmp <- spe[x, ]
  tmp <- tmp[, colSums(logcounts(tmp)) > 0]
  imgData(tmp) <- NULL
  tmp
})
cat("\nAdditional sampling with recalculated n=5000 logcounts to make consistent with MBv implementation\n")
tmp <- computeLibraryFactors(speList[["n5000"]])
tmp <- logNormCounts(tmp)
tmp <- tmp[, colSums(logcounts(tmp)) > 0]
speList$n5000_relognorm <- tmp
cat("\nsave speList",format(Sys.time(), tz="EST"),"\n")
sapply(speList, dim)
saveRDS(speList, "processed-data/04_feature_selection/spoon-tutorial_test/speList.rda")

cat("\nnnSVG",format(Sys.time(), tz="EST"),"\n")
nnSVG_results <- lapply(names(speList), function(x) {
  set.seed(123)
  rowData(nnSVG(speList[[x]], assay="logcounts", n_threads=12))
})
names(nnSVG_results) <- names(speList)
cat("\nsave nnSVG_results",format(Sys.time(), tz="EST"),"\n")
saveRDS(nnSVG_results, "processed-data/04_feature_selection/spoon-tutorial_test/nnSVG_results.rda")

cat("\nweightList",format(Sys.time(), tz="EST"),"\n")
weightList <- lapply(names(speList), function(x) {
  set.seed(123)
  generate_weights(input = speList[[x]], stabilize = TRUE, n_threads=12)
})
names(weightList) <- names(speList)
cat("\nsave weightList",format(Sys.time(), tz="EST"),"\n")
saveRDS(weightList, "processed-data/04_feature_selection/spoon-tutorial_test/weightList.rda")

source("code/04_feature_selection/spoon-tutorial_test/weighted-nnSVG_re-write_iterative-cov-matrix.r")
cat("\nweightedJT nnSVG",format(Sys.time(), tz="EST"),"\n")
w.nnSVG_results <- lapply(names(speList), function(x) {
  tmp = speList[[x]]
  weighted_logcounts <- t(weightList[[x]])*assays(tmp)[['logcounts']]
  assay(tmp, "weighted_logcounts") <- weighted_logcounts
  set.seed(123)
  rowData(weightedJT_nnSVG(tmp, X=weightList[[x]], assay_name="weighted_logcounts", n_threads=12))
})
names(w.nnSVG_results) <- names(speList)
cat("\nsave w.nnSVG_results",format(Sys.time(), tz="EST"),"\n")
saveRDS(w.nnSVG_results, "processed-data/04_feature_selection/spoon-tutorial_test/w.nnSVG_results.rda")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
