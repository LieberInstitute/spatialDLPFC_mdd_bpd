setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(scran)
	library(here)
})

load(here("processed-data","03_QC","spe_demo-filt.Rdata"))

spe$exclude_outliers = spe$umi_local.outlier | spe$genes_local.outlier | spe$chrM.ratio_local.outlier | spe$exclude_east_edge
spe$extra_3MAD_outliers = spe$umi_3MAD.outlier_sample | spe$genes_3MAD.outlier_sample
spe$plot_outliers = factor(paste(spe$exclude_outliers, spe$extra_3MAD_outliers),
	levels=c("FALSE FALSE","FALSE TRUE","TRUE FALSE","TRUE TRUE"),
	labels=c("keep","3MAD flag","exclude","exclude"))
(spe = spe[,spe$plot_outliers!="exclude"])

spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)

save(spe, file=here("processed-data","03_QC","spe_demo-filt.Rdata"))

write(c(paste("Filtered spe_demo-filt on",format(Sys.time(), tz="UTC"),"UTC"),
	paste("File location:",here("processed-data","03_QC","spe_demo-filt.Rdata")),
	paste("Source code:",here("code","03_QC","04_filter-outliers.r")),
	"*","*","*"), here("spe_tracker_current.txt"), append=TRUE)
