setwd("/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd")
library(SpatialExperiment)
library(scran)
library(STexampleData)
library(spoon)
library(BiocParallel)
library(nnSVG)

spe <- Visium_mouseCoronal()
spe <- spe[, colData(spe)$in_tissue == 1]
spe <-nnSVG:: filter_genes(spe)
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)

set.seed(123)
ix_random <- sample(seq_len(nrow(spe)), 10)
spe <- spe[ix_random, ]
spe <- spe[, colSums(logcounts(spe)) > 0]
table(rowSums(logcounts(spe))==0)
table(colSums(logcounts(spe))==0)

table(rowSums(counts(spe))==0)
table(colSums(counts(spe))==0)

#test reproducibility of nnSVG
spe1 = nnSVG(spe, assay="logcounts")
spe2 = nnSVG(spe, assay="logcounts")

###no
plot(rowData(spe1)$tau.sq, rowData(spe2)$tau.sq)
abline(a=0, b=1)

#test reproducibility with re-setting global seed before hand 
set.seed(123)
spe3 = nnSVG(spe, assay="logcounts")
set.seed(123)
spe4 = nnSVG(spe, assay="logcounts")

###yes
plot(rowData(spe3)$tau.sq, rowData(spe4)$tau.sq)
abline(a=0, b=1)

#test if using BPPARAM seed is also reproducible
spe5 = nnSVG(spe, assay="logcounts", BPPARAM=MulticoreParam(RNGseed=12))
spe6 = nnSVG(spe, assay="logcounts", BPPARAM=MulticoreParam(RNGseed=12))

###no (meaning the parallelization of BRISC_estimation (the only place where bplapply occurs and BPPARAM is passed), is not influenced by RNGseed)
plot(rowData(spe5)$tau.sq, rowData(spe6)$tau.sq)
abline(a=0, b=1)




#calculate weights for weighted nnSVG tests
set.seed(123)
weights <- generate_weights(input = spe, stabilize = TRUE)

#test reproducibility of my modified nnSVG
source("code/04_feature_selection/spoon-tutorial_test/weighted-nnSVG_re-write_iterative-cov-matrix.r")
weighted_logcounts <- t(weights)*assays(spe)[['logcounts']]
weighted_mean <- Matrix::rowMeans(weighted_logcounts)
spe_JT = spe
assay(spe_JT, "weighted_logcounts") <- weighted_logcounts

set.seed(123)
spe_JT1 <- weightedJT_nnSVG(spe_JT, X=weights, assay_name="weighted_logcounts")
set.seed(123)
spe_JT2 <- weightedJT_nnSVG(spe_JT, X=weights, assay_name="weighted_logcounts")

###yes
plot(rowData(spe_JT1)$weighted_tau.sq, rowData(spe_JT2)$weighted_tau.sq)
abline(a=0, b=1)


#test reproducibility of KS weighted_nnSVG
set.seed(123)
spe_KS1 <- weighted_nnSVG(input = spe, w = weights)
#get low/zero expr gene error
set.seed(123)
spe_KS2 <- weighted_nnSVG(input = spe, w = weights)
#get low/zero expr gene error

###no
plot(rowData(spe_KS1)$weighted_tau.sq, rowData(spe_KS2)$weighted_tau.sq)
abline(a=0, b=1)

#test reproducibility with BPPARAM seed only
spe_KS3 <- weighted_nnSVG(input = spe, w = weights, BPPARAM=MulticoreParam(RNGseed = 12))
spe_KS4 <- weighted_nnSVG(input = spe, w = weights, BPPARAM=MulticoreParam(RNGseed = 12))

###yes
plot(rowData(spe_KS3)$weighted_tau.sq, rowData(spe_KS4)$weighted_tau.sq)
abline(a=0, b=1)

#so what if I set RNGseed to be 123, will that match the results from my weighted nnSVG?
spe_KS5 <- weighted_nnSVG(input = spe, w = weights, BPPARAM=MulticoreParam(RNGseed = 123))
spe_KS6 <- weighted_nnSVG(input = spe, w = weights, BPPARAM=MulticoreParam(RNGseed = 123))

###match each other
plot(rowData(spe_KS5)$weighted_tau.sq, rowData(spe_KS6)$weighted_tau.sq)
abline(a=0, b=1)

###match me? -- no
plot(rowData(spe_KS5)$weighted_tau.sq, rowData(spe_JT1)$weighted_tau.sq)
abline(a=0, b=1)

###match me mean? -- yes (sanity check)
plot(rowData(spe_KS5)$weighted_mean, rowData(spe_JT1)$weighted_mean)
abline(a=0, b=1)

#set all seeds for both?
set.seed(123)
spe_JT3 <- weightedJT_nnSVG(spe_JT, X=weights, assay_name="weighted_logcounts", BPPARAM=MulticoreParam(RNGseed = 12))
set.seed(123)
spe_KS7 <- weighted_nnSVG(input = spe, w = weights, BPPARAM=MulticoreParam(RNGseed = 12))

###match me? -- no
plot(rowData(spe_KS7)$weighted_tau.sq, rowData(spe_JT3)$weighted_tau.sq)
abline(a=0, b=1)

###match me mean? -- yes (sanity check)
plot(rowData(spe_KS7)$weighted_mean, rowData(spe_JT3)$weighted_mean)
abline(a=0, b=1)

###match prev run? -- yes for both (because RNGseed is disposable for BRISC_estimation in my function)
plot(rowData(spe_KS3)$weighted_tau.sq, rowData(spe_KS7)$weighted_tau.sq)
abline(a=0, b=1)

plot(rowData(spe_JT1)$weighted_tau.sq, rowData(spe_JT3)$weighted_tau.sq)
abline(a=0, b=1)


###match prev run with different RNGseed -- no (good, confirmatory)
plot(rowData(spe_KS5)$weighted_tau.sq, rowData(spe_KS7)$weighted_tau.sq)
abline(a=0, b=1)


#compare magnitude of difference between seeds to that between KS and JT method
### tau.sq
plot(rowData(spe_KS5)$weighted_tau.sq, rowData(spe_KS7)$weighted_tau.sq)
abline(a=0, b=1)
plot(rowData(spe_JT3)$weighted_tau.sq, rowData(spe_KS7)$weighted_tau.sq)
abline(a=0, b=1)

### sigma.sq
plot(rowData(spe_KS5)$weighted_sigma.sq, rowData(spe_KS7)$weighted_sigma.sq)
abline(a=0, b=1)
plot(rowData(spe_JT3)$weighted_sigma.sq, rowData(spe_KS7)$weighted_sigma.sq)
abline(a=0, b=1)

### LR stat
plot(rowData(spe_KS5)$weighted_LR_stat, rowData(spe_KS7)$weighted_LR_stat)
abline(a=0, b=1)
plot(rowData(spe_JT3)$weighted_LR_stat, rowData(spe_KS7)$weighted_LR_stat)
abline(a=0, b=1)


# CONCLUSION: without modifying KS/spoon code, we won't be able to get the exact same results due to seed/ differences in looping internal BRISC_estimation vs looping nnSVG function itself
### however, the differences between my modified code and KS code are on the same scale as the differences in tau.sq due to seeds, meaning that the two functions are statistically equivalent

# ADDITIONAL NOTE: KS correctly modified the weighted_mean calc to occur outside of nnSVG because nnSVG internally forces the mean and a few other stats to be calculated on the logcounts assay
### she didn't do this for the other stat, so var is not accurate
### https://github.com/lmweber/nnSVG/blob/devel/R/nnSVG.R#L282 ### code line where nnSVG forces mean and var calc from logcounts not the named assay

