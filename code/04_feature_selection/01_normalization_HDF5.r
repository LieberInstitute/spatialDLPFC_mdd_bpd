setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(scuttle)
})
set.seed(123)
setAutoBlockSize(1e9)

cdata = read.csv("processed-data/03_QC/colData_edges-problem-areas_spotsweeper_FINAL.csv", row.names=1)

system.time(spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_"))
cat("\nDim spe:",dim(spe),"\n")
#cat("\ntree:\n")
#showtree(counts(spe))
#cat("\nseed:\n")
#seed(counts(spe))

stopifnot(identical(rownames(cdata), rownames(colData(spe))))

spe$remove_spots = cdata$remove_spots
spe$problem_area_flag = cdata$problem_area_flag

format(Sys.time(), tz="EST")
spe = spe[,spe$remove_spots==FALSE]
cat("\nDim spe (QC filtered):",dim(spe),"\n\n")

spe = spe[,spe$sample_id!="V13F27-338_C1"]
cat("\nDim spe (extra Br5366 sample removed):",dim(spe),"\n\n")

format(Sys.time(), tz="EST")
spots_more0 = rowSums(counts(spe)>0)
spe = spe[spots_more0>5,]
cat("\nDim spe (genes with non-zero count in <=5 spots removed):",dim(spe),"\n\n")

#rearrange
fix.order = order(paste(spe$round, spe$sample_id, spe$array_col, spe$array_row))
spe = spe[,fix.order]

#normalization
format(Sys.time(), tz="EST")
cat("\nCompute library factors and normalize counts")
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)

#save HDF5
start.time = Sys.time()
cat("\nStart HDF5SummarizedExperiment save:"); start.time

saveHDF5SummarizedExperiment(spe, dir="processed-data/04_feature_selection", prefix="spe_n120_postQC_norm_",
        chunkdim=c(100,500),
        verbose=F)

cat("\nTime elapsed (saveHDF5):",
        round(difftime(Sys.time(), start.time, units="hours"),2), "hours\n")

#update spe tracker
write(c(paste("**** Created filtered, normalized spe on",format(Sys.time(), tz="EST"),"EST"),
        "**** Old file location: processed-data/02_build_spe/spe_n120_",
        "**** New file location: processed-data/04_feature_selection/spe_n120_postQC_norm_",
        "**** Source code: code/04_feature_selection/01_normalization_HDF5.r",
        "****","****","****"), "spe_tracker_current.txt", append=TRUE)

## Reproducibility information
cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
