suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
})
set.seed(123)

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n120_postQC_norm_")
cat("Dim spe\n")
dim(spe)

#orig chunk dim is 1 x n obs
cat("\n\nauto block size:\n")
getAutoBlockSize()
cat("\n\noriginal chunk dim:\n")
chunkdim(counts(spe))
cat("\nrowMax implementation",format(Sys.time(), tz="EST"),"...\n")
system.time(rowMaxs(counts(spe)))

#increase block size 
cat("\n\nincrease block size:\n")
setAutoBlockSize(1e9)
cat("\n\noriginal chunk dim:\n")
chunkdim(counts(spe))
cat("\nrowMax implementation",format(Sys.time(), tz="EST"),"...\n")
system.time(rowMaxs(counts(spe)))

#default chunk dim
cat("\n\ntry default chunk dim:\n")
dim(defaultAutoGrid(counts(spe)))
cat("\nRealize new chunk dims",format(Sys.time(), tz="EST"),"...\n")
dim2 = writeHDF5Array(counts(spe), chunkdim=dim(defaultAutoGrid(counts(spe))))
cat("rowMax implementation",format(Sys.time(), tz="EST"),"...\n")
system.time(rowMaxs(dim2))

#set chunk dim
cat("\n\ntry manual chunk dim:\n")
c(100,500)
cat("\nRealize new chunk dims",format(Sys.time(), tz="EST"),"...\n")
dim3 = writeHDF5Array(counts(spe), chunkdim=c(100,500))
cat("rowMax implementation",format(Sys.time(), tz="EST"),"...\n")
system.time(rowMaxs(dim3))

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
