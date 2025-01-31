setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(spdep)
	library(scater)
	library(SpotSweeper)
	library(here)
})
set.seed(123)

#load spe
start.time = Sys.time()
cat("Start time (loadHDF5SummarizedExperiment):", format(start.time),"\n")
spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_")
cat("Time elapsed (loadHDF5SummarizedExperiment):", round(difftime(Sys.time(), start.time, units="mins"),2), "minutes\n")
#load coldata
cdata = read.csv(here("processed-data","03_QC","spe_n120_edge-detection_colData.csv"), row.names=1)
stopifnot(identical(rownames(colData(spe)), rownames(cdata)))
#add edge info and subset
colData(spe)$keep_spots = cdata$keep_spots
colData(spe)$edge_outlier = cdata$edge_outlier
cat("dim spe before edge outlier removal:",dim(spe),"\n")
spe = spe[,spe$edge_outlier==FALSE & spe$in_tissue==TRUE]
cat("dim spe after edge outlier removal:",dim(spe),"\n\n")

#prep for outlier detection
colData(spe)$lg10.umi = log10(colData(spe)$sum_umi)
colData(spe)$lg10.genes = log10(colData(spe)$sum_gene)

#run spotsweeper
cat(format(Sys.time()), "Calculate local outliers (umi counts)...","\n")
spe <- localOutliers(spe, metric = "sum_umi", direction = "lower", log = TRUE)
#colData(spe)[,c(ncol(colData(spe)),ncol(colData(spe))-2)] = NULL
colnames(colData(spe))[ncol(colData(spe))-1] = "umi_local.outlier"

cat(format(Sys.time()), "Calculate local outliers (n genes)...","\n")
spe <- localOutliers(spe, metric = "sum_gene", direction = "lower", log = TRUE)
#colData(spe)[,c(ncol(colData(spe)),ncol(colData(spe))-2)] = NULL
colnames(colData(spe))[ncol(colData(spe))-1] = "genes_local.outlier"

cat(format(Sys.time()), "Calculate local outliers (mito %)...","\n")
spe <- localOutliers(spe, metric = "expr_chrM_ratio", direction = "higher", log = FALSE)
#colData(spe)[,c(ncol(colData(spe)),ncol(colData(spe))-2)] = NULL
colnames(colData(spe))[ncol(colData(spe))-1] = "chrM.ratio_local.outlier"

#update coldata
true.spots = cdata$edge_outlier==FALSE & cdata$in_tissue==TRUE
stopifnot(identical(rownames(cdata)[true.spots], rownames(colData(spe))))
cdata$umi_local.outlier = NA
cdata[true.spots, "umi_local.outlier"] = colData(spe)$umi_local.outlier
cdata$genes_local.outlier = NA
cdata[true.spots, "genes_local.outlier"] = colData(spe)$genes_local.outlier
cdata$chrM.ratio_local.outlier = NA
cdata[true.spots, "chrM.ratio_local.outlier"] =	colData(spe)$chrM.ratio_local.outlier

#save coldata
write.csv(cdata, here("processed-data","03_QC","spe_n120_edge-detection_spotsweeper_colData.csv"), row.names=T)
cat("\n", format(Sys.time()), "Updated colData saved to:", here("processed-data","03_QC","spe_n120_edge-detection_spotsweeper_colData.csv"),"\n")
#save(spe, file=here("processed-data","03_QC","spe_demo-filt.Rdata"))

#write(c(paste("******* Modified spe_demo-filt on",format(Sys.time(), tz="UTC"),"UTC"),
#        paste("******* Old location:",here("processed-data","03_QC","spe_demo-filt.Rdata")),
#	paste("******* New location:",here("processed-data","03_QC","spe_demo-filt.Rdata")),
#        paste("******* Source code:",here("code","03_QC","02_outliers.r")),
#        "*********","*********","*********"), here("spe_tracker_current.txt"), append=TRUE)
#
## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
