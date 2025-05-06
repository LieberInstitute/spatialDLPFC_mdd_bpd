setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(edgeR)
})
set.seed(123)
#setAutoBlockSize(1e9)

#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
load("processed-data/06_pseudobulk/spe_n119_pseudo_sample-comb-clus_norm-filt.Rdata")
spe_pseudo$slide_name = gsub("-","\\.", spe_pseudo$slide)
#spe_pseudo$lg10_sum = log10(spe_pseudo$sum) #can't include sum as covar without log10 transformation without breaking fullrank==T

### the way that makes more sense to me for setting up contrasts
dx = spe_pseudo$condition
clus = spe_pseudo$combined_cluster
sex = spe_pseudo$sex
group = interaction(dx, clus, sex)
#group = interaction(dx, clus)

dx_mod <- model.matrix(
  ~ 0 + group + age + slide_name + detected,
  colData(spe_pseudo)
)

#dx_mod <- model.matrix(
#  ~ 0 + group + age + sex + slide_name,
#  colData(spe_pseudo)
#)
stopifnot(is.fullrank(dx_mod))

corfit <- duplicateCorrelation(
  logcounts(spe_pseudo),
  design = dx_mod,
  block = colData(spe_pseudo)$sample_id
)

fit <- lmFit(
  logcounts(spe_pseudo),
  design = dx_mod,
  block = colData(spe_pseudo)$sample_id,
  correlation = corfit$consensus
)

#colnames(coef(fit))
#sexes combined, by condition within sex
str1 = as.character(outer(levels(spe_pseudo$combined_cluster), c("NTC.MDD","NTC.BPD","MDD.BPD"), paste, sep="_"))
str2 = as.character(outer(str1, c("both","M","F"), paste, sep="_"))
#by sex within condition
str3 = as.character(outer(levels(spe_pseudo$combined_cluster), c("NTC.NTC_sex","MDD.MDD_sex","BPD.BPD_sex"), paste, sep="_"))

comparisons = c(str2, str3)
cont_mtx.s = matrix(0, nrow=ncol(coef(fit)), ncol=length(comparisons), dimnames = list(colnames(coef(fit)),comparisons))

#comparisons = as.character(outer(levels(spe_pseudo$combined_cluster), c("NTC.MDD_none","NTC.BPD_none","MDD.BPD_none"), paste, sep="_"))
#cont_mtx = matrix(0, nrow=ncol(coef(fit)), ncol=length(comparisons), dimnames = list(colnames(coef(fit)),comparisons))

for(i in comparisons) {
  seg = unlist(strsplit(i,"_"))
  clus = seg[[1]]
  dx_refer = unlist(strsplit(seg[[2]], "\\."))[[1]]
  dx_compare = unlist(strsplit(seg[[2]], "\\."))[[2]]
#  if(seg[[3]]=="none") {
#    cont_mtx[grep(paste(dx_refer, clus, sep="."), rownames(cont_mtx)),i] = -1
#    cont_mtx[grep(paste(dx_compare, clus, sep="."), rownames(cont_mtx)),i] = 1
#  }
  if(seg[[3]]=="both") {
    cont_mtx.s[grep(paste(dx_refer, clus,sep="."), rownames(cont_mtx.s)),i] = c(-1,-1)
    cont_mtx.s[grep(paste(dx_compare, clus,sep="."), rownames(cont_mtx.s)),i] = c(1,1)
  }
  if(seg[[3]]=="sex") {
    cont_mtx.s[grep(paste(dx_refer, clus, "M", sep="."), rownames(cont_mtx.s)),i] = -1
    cont_mtx.s[grep(paste(dx_compare, clus, "F", sep="."), rownames(cont_mtx.s)),i] = 1
  }
  if(seg[[3]]=="M") {
    cont_mtx.s[grep(paste(dx_refer, clus, "M", sep="."), rownames(cont_mtx.s)),i] = -1
    cont_mtx.s[grep(paste(dx_compare, clus, "M", sep="."), rownames(cont_mtx.s)),i] = 1
  }
  if(seg[[3]]=="F") {
    cont_mtx.s[grep(paste(dx_refer, clus, "F", sep="."), rownames(cont_mtx.s)),i] = -1
    cont_mtx.s[grep(paste(dx_compare, clus, "F", sep="."), rownames(cont_mtx.s)),i] = 1
  }
}

tmp1 = contrasts.fit(fit, cont_mtx.s)
#tmp1 = contrasts.fit(fit, cont_mtx)
tmp2 = eBayes(tmp1)

saveRDS(tmp2, "processed-data/07_dx_DE/eBayes_fit_layer-specific_condition-sex_covars-age-detected-slide.rda")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
