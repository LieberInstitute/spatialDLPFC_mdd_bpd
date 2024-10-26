setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
        library(SpatialExperiment)
        library(HDF5Array)
        library(scuttle)
        library(here)
})
set.seed(123)

#load spe
start.time = Sys.time()
cat("Start time (loadHDF5SummarizedExperiment):", format(start.time),"\n")
spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_")
cat("Time elapsed (loadHDF5SummarizedExperiment):", round(difftime(Sys.time(), start.time, units="mins"),2), "minutes\n")

#load column data
cdata = read.csv(here("processed-data","03_QC","spe_n120_edge-detection_spotsweeper_colData.csv"), row.names=1)
stopifnot(identical(rownames(colData(spe)), rownames(cdata)))

#add edge info and subset
cat("\n# spots off tissue:\n")
table(!spe$in_tissue)
colData(spe)$keep_spots = cdata$keep_spots
colData(spe)$edge_outlier = cdata$edge_outlier
cat("\n# spots edge outliers:\n")
table(spe$edge_outlier)
colData(spe)$local_outlier = cdata$umi_local.outlier | cdata$genes_local.outlier | cdata$chrM.ratio_local.outlier
cat("\n# spots local outliers:\n")
table(spe$local_outlier)

#remove outliers
any_outlier = spe$edge_outlier | spe$local_outlier | !spe$in_tissue
table(any_outlier)
cat("dim spe before outlier removal:",dim(spe),"\n")
spe = spe[,any_outlier==FALSE]
cat("dim spe after outlier removal:",dim(spe),"\n\n")
spe = spe[rowSums(counts(spe))!=0,]
cat("dim spe after zero-count genes removal:",dim(spe),"\n\n")

#normalization
#spe <- computeLibraryFactors(spe)
#spe <- logNormCounts(spe)

#save HDF5
start.time2 = Sys.time()
cat("\nStart save:"); start.time2

saveHDF5SummarizedExperiment(spe, dir=here("processed-data","03_QC"), prefix="spe_n120_postQC_",
        chunkdim=getHDF5DumpChunkDim(c(1,ncol(spe))),
        verbose=T)

cat("\nTime elapsed (saveHDF5):",
        round(difftime(Sys.time(), start.time2, units="hours"),2), "hours\n")

#update spe tracker
write(c(paste("*** Created post QC filtered spe on",format(Sys.time(), tz="UTC"),"UTC"),
        paste("*** Old file location:",here("processed-data","02_build_spe","spe_n120_")),
        paste("*** New file location:",here("processed-data","03_QC","spe_n120_postQC_")),
        paste("**** Source code:",here("code","03_QC","03_filter-outliers_HDF5.r")),
        "***","***","***"), here("spe_tracker_current.txt"), append=TRUE)

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
