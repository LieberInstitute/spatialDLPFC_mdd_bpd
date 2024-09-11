setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
    library("here")
    library("SpatialExperiment")
    library("spatialLIBD")
    library("rtracklayer")
    library("HDF5Array")
})

## Define some info for the samples
fastqs = grep("^V", list.files('raw-data/FASTQ'), value=T)
space = grep("^V", list.files("processed-data/01_spaceranger"), value=T)

stopifnot(identical(fastqs, space))

proc_slide = sapply(strsplit(space,"_"), function(x) x[1])
proc_array = sapply(strsplit(space,"_"), function(x) x[2])

#fix mbv demo table
demo = read.csv("raw-data/sample_info/MBv_demographics_combined.csv")
### first change 8 to B
demo$slide = gsub("V138","V13B",demo$slide)
### then change the 1 sample that had a different image number
demo[demo$mbv_sample=="MBv_081","slide"] = "V13B23-283"
demo[demo$mbv_sample=="MBv_081","array"] = "B1"
#then assign the samples missing slide info to the remaining slide
missing_slide = setdiff(proc_slide, demo$slide)
demo[is.na(demo$slide),"slide"] = missing_slide

#make sure that all samples with output data are represented in the demo table
demo$sample_id = paste(demo$slide, demo$array, sep="_")
stopifnot(identical(sort(space),sort(demo$sample_id)))

demo$sample_path = file.path(here("processed-data", "01_spaceranger"), demo$sample_id,"outs")
#sample_paths = file.path(here::here("processed-data", "01_spaceranger"), demo$sample_id,"outs")

## Build basic SPE
start.time = Sys.time()
cat("Start time:"); start.time; cat("\n")
spe <- read10xVisiumWrapper(
	demo$sample_path,
	demo$sample_id,
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
	new_col <- merge(colData(spe), demo)
	## Fix order
	new_col <- new_col[match(spe$key, new_col$key), ]
	stopifnot(identical(new_col$key, spe$key))
	rownames(new_col) <- rownames(colData(spe))
	colData(spe) <- new_col[, -which(colnames(new_col) == "sample_path")]
	return(spe)
}
spe <- add_design(spe)

#change colnames to be completely unique
spotcodes = paste(rownames(colData(spe)), spe$sample_id, sep="_")
colnames(spe) <- spotcodes
rownames(colData(spe)) <- spotcodes

#save
start.time2 = Sys.time()
cat("\nStart save:"); start.time2           

saveHDF5SummarizedExperiment(spe, dir=here("processed-data","02_build_spe"), prefix="spe_n120_",
        chunkdim=getHDF5DumpChunkDim(c(1,ncol(spe))),
        verbose=F)

cat("\nTime elapsed (saveHDF5):",
        round(difftime(Sys.time(), start.time2, units="mins"),2), "minutes\n")

if(file.exists(here("spe_tracker_current.txt"))) {
	if(file.exists(here("spe_tracker_archive.txt"))) {
		x = readLines(here("spe_tracker_current.txt"))
		write(c("########", paste("ARCHIVED",format(Sys.time(), tz="UTC"),"UTC"), "########"), here("spe_tracker_archive.txt"), append=TRUE)
	} else {
		writeLines(c("########", paste("ARCHIVED",format(Sys.time(), tz="UTC"),"UTC"), "########"), file(here("spe_tracker_archive.txt"))) 
		x = readLines(here("spe_tracker_current.txt"))
		write(x, here("spe_tracker_archive.txt"), append=TRUE)
		close(file(here("spe_tracker_archive.txt")))
	}
}
writeLines(c(paste("Created raw spe (HDF5) on",format(Sys.time(), tz="UTC"),"UTC"), 
	paste("New file location:",here("processed-data","02_build_spe","spe_n120_")), 
	paste("Source code:",here("code","02_build_spe","01_raw_spe_hdf5.r")),"*","*","*"), 
	file(here("spe_tracker_current.txt")))
close(file(here("spe_tracker_current.txt")))

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
