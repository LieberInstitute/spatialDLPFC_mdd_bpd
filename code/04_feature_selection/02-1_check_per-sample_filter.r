setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(parallel)
})
set.seed(123)

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n120_postQC_norm_")
#nnSVG vignette says to filter each sample individually: https://bioconductor.uib.no/packages/3.20/bioc/vignettes/nnSVG/inst/doc/nnSVG.Rmd
### when i was running nnSVG per slide, I had to run 4 samples/each slide on the same set of genes
### so minimally I had to filter each slide individually or collapse the per-sample filters up to a per-slide level
### since this had to happen anyways, i thought it would be better to run all slides on the same group of genes 
### however, now that i am running each sample individually, some samples are erroring because there are too many high 0 genes 
### sooooo it is probably better to filter each sample independently
### but how does that impact the genes that get included?

sample.list = unique(spe$sample_id)
names(sample.list) = sample.list

setAutoBlockSize(1e9)

#alt1 = DelayedArray::rowMaxs(counts(spe)) 
alt1 = DelayedArray::rowMaxs(logcounts(spe))
table(alt1>2) #can first remove genes that don't even have a max count of >2 in any spot in any sample

keep.genes = names(alt1)[alt1>4] #increased to >4 for logcounts after >2 was too loose
length(keep.genes)
#t1 = counts(spe)[keep.genes,]
t1 = logcounts(spe)[keep.genes,]
dim(t1)

l1 = mclapply(sample.list, function(x) {
  t2 = t1[,spe$sample_id==x]
  t3 = rowSums(t2>4) #increased to >4 for logcounts after >2 was too loose
  #t4 = t3>=5
  t4 = t3>=10 #increased to >=10 for logcounts
  names(t4)[t4]
}, mc.cores = 12 
)

sapply(l1, length)
saveRDS(l1, "processed-data/04_feature_selection/per-sample_gene-filter_logcounts.rda")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
