setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
        library(SpatialExperiment)
        library(here)
})

load(here("processed-data","03_QC","spe_demo-filt.Rdata"))

exclude_outliers = spe$umi_local.outlier | spe$genes_local.outlier | spe$chrM.ratio_local.outlier
cat("Pre-filter dimensions:",dim(spe),"\n")
(spe = spe[,!exclude_outliers])
cat("Dim. after local outlier filter:",dim(spe),"\n")

save(spe, file=here("processed-data","03_QC","spe_demo-filt-outliers.Rdata"))

write(c(paste("********* Filtered spe_demo-filt on",format(Sys.time(), tz="UTC"),"UTC"),
        paste("********* Old file location:",here("processed-data","03_QC","spe_demo-filt.Rdata")),
        paste("********* New file location:",here("processed-data","03_QC","spe_demo-filt-outliers.Rdata")),
	paste("********* Source code:",here("code","03_QC","03_filter-outliers.r")),
        "***********","***********","***********"), here("spe_tracker_current.txt"), append=TRUE)
