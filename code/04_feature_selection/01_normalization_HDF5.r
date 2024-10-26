setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(scuttle)
})
set.seed(123)
#DelayedArray:::set_verbose_block_processing(TRUE)

spe <- loadHDF5SummarizedExperiment(dir="processed-data/03_QC/", prefix="spe_n120_postQC_")
cat("\n\nspe dimensions:\n")
dim(spe)
cat("\ntree:\n")
showtree(counts(spe))
cat("\nseed:\n")
seed(counts(spe))

cat("\nnset verbose block processing:\n")
DelayedArray:::set_verbose_block_processing(TRUE)

cat("\n\nCompute library factors and normalize counts")
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)


#save backups
cat("\nwriteHDF5Array counts and logcounts backup\n")
cmtx = as(counts(spe), "HDF5Matrix")
writeHDF5Array(cmtx, filepath="processed-data/04_feature_selection/spe_n120_postQC_raw-counts_backup",
               name="counts")

lcmtx = as(logcounts(spe), "HDF5Matrix")
writeHDF5Array(lcmtx, filepath="processed-data/04_feature_selection/spe_n120_postQC_logcounts_backup",
               name="logcounts")

#save HDF5
start.time = Sys.time()
cat("\nStart HDF5SummarizedExperiment save:"); start.time

saveHDF5SummarizedExperiment(spe, dir=here("processed-data","04_feature_selection"), prefix="spe_n120_postQC_norm_",
        chunkdim=getHDF5DumpChunkDim(c(1,ncol(spe))),
        verbose=T)

cat("\nTime elapsed (saveHDF5):",
        round(difftime(Sys.time(), start.time, units="hours"),2), "hours\n")

#update spe tracker
write(c(paste("****** Created normalized spe on",format(Sys.time(), tz="UTC"),"UTC"),
        paste("****** Old file location:",here("processed-data","03_QC","spe_n120_postQC_")),
        paste("****** New file location:",here("processed-data","04_feature_selection","spe_n120_postQC_norm_")),
        paste("****** Source code:",here("code","04_feature_selection","01_normalization_HDF5.r")),
        "******","******","******"), here("spe_tracker_current.txt"), append=TRUE)

## Reproducibility information
cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
