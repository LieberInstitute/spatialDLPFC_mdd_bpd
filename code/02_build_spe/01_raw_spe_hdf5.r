
setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
    library("here")
    library("SpatialExperiment")
    library("spatialLIBD")
    library("rtracklayer")
    library("HDF5Array")
})

## Define some info for the samples
load(here::here("code", "REDCap", "REDCap_MBv.rda"))

sample_info <- data.frame(slide = as.factor(REDCap_MBv$slide))
sample_info$array <- as.factor(REDCap_MBv$array)
sample_info$brnum <- as.factor(sapply(strsplit(REDCap_MBv$sample, "-"), `[`, 1))
sample_info$species <- as.factor(REDCap_MBv$species)
sample_info$replicate <- as.factor(REDCap_MBv$serial)
sample_info$sample_id <- paste(sample_info$slide, sample_info$array, sep = "_")
sample_info$sample_path = file.path(here::here("processed-data", "01_spaceranger"), sample_info$sample_id,"outs")

list4spe = c('V13F27-338','V13F27-348', 'V13Y10-020','V13Y10-021', 'V13Y10-022','V13Y10-023')
sample_info = sample_info[sample_info$slide %in% list4spe,]

## Build basic SPE
start.time = Sys.time()
cat("Start time:"); start.time; cat("\n")
spe <- read10xVisiumWrapper(
    sample_info$sample_path,
    sample_info$sample_id,
    type = "HDF5",
    data = "raw",
    images = c("lowres"),
    load = FALSE,
    reference_gtf = file.path("/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/","genes", "genes.gtf")
)
cat("\nTime elapsed (read10x, HDF5):",
	round(difftime(Sys.time(), start.time, units="mins"),2), "minutes\n")

start.time2 = Sys.time()
cat("\nStart save:"); start.time2

## Add the study design info
add_design <- function(spe) {
    new_col <- merge(colData(spe), sample_info)
    ## Fix order
    new_col <- new_col[match(spe$key, new_col$key), ]
    stopifnot(identical(new_col$key, spe$key))
    rownames(new_col) <- rownames(colData(spe))
    colData(spe) <-
        new_col[, -which(colnames(new_col) == "sample_path")]
    return(spe)
}
spe <- add_design(spe)

#change colnames to be completely unique
spotcodes = paste(spe$sample_id, rownames(colData(spe)), "_")
colnames(spe) <- spotcodes
rownames(colData(spe)) <- spotcodes

#save
start.time2 = Sys.time()
cat("\nStart save:"); start.time2           

saveHDF5SummarizedExperiment(spe, dir=here("processed-data","02_build_spe"), prefix="test_spe_2d_",
        chunkdim=getHDF5DumpChunkDim(c(1,ncol(spe))),
        verbose=F)

cat("\nTime elapsed (saveHDF5):",
        round(difftime(Sys.time(), start.time2, units="mins"),2), "minutes\n")

#if(file.exists(here("spe_tracker_current.txt"))) {
#	if(file.exists(here("spe_tracker_archive.txt"))) {
#		x = readLines(here("spe_tracker_current.txt"))
#		write(c("########", paste("ARCHIVED",format(Sys.time(), tz="UTC"),"UTC"), "########"), here("spe_tracker_archive.txt"), append=TRUE)
#	} else {
#		writeLines(c("########", paste("ARCHIVED",format(Sys.time(), tz="UTC"),"UTC"), "########"), file(here("spe_tracker_archive.txt"))) 
#		x = readLines(here("spe_tracker_current.txt"))
#		write(x, here("spe_tracker_archive.txt"), append=TRUE)
#		close(file(here("spe_tracker_archive.txt")))
#	}
#}
#writeLines(c(paste("Created spe_raw on",format(Sys.time(), tz="UTC"),"UTC"), 
#	paste("New file location:",here("processed-data","02_build_spe","spe_raw.Rdata")), 
#	paste("Source code:",here("code","02_build_spe","01_raw_spe.R")),"*","*","*"), 
#	file(here("spe_tracker_current.txt")))
#close(file(here("spe_tracker_current.txt")))

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
