setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
    library("SpatialExperiment")
    library("spatialLIBD")
    library("rtracklayer")
    library("HDF5Array")
})
source("code/02_build_spe/getMBvSampleInfo_function.r")

## Identify samples
fastqs = grep("^V", list.files('raw-data/FASTQ'), value=T)
space = grep("^V", list.files("processed-data/01_spaceranger"), value=T)

stopifnot(identical(fastqs, space))

## Get sample info
demo = getMBvSampleInfo(REDCapFile="Visium_DATA_2025-01-22_1406.csv",
                        demoFile="DLPFC_cross-disorders_demographics_MBv.csv")

#list any samples from demo data that are missing spaceranger output
missingSpaceranger = setdiff(demo$sample_id, space)
cat("\nSamples from demo/exp data that are missing spaceranger output.\n(should only be V13B23-339_A1)\n")
missingSpaceranger

stopifnot(missingSpaceranger=="V13B23-339_A1")
demo = demo[demo$sample_id!="V13B23-339_A1",]

#make sure that all samples with output data are represented in the demo table
stopifnot(identical(sort(space),sort(demo$sample_id)))
demo$sample_path = file.path("processed-data/01_spaceranger", demo$sample_id,"outs")

#now remove unwanted samples/slides from the list used to generate spe
cat("\nRemove samples from slides V13Y10-020 and V13B23-331 from the list used to generate the spe\n")
demo = demo[!demo$slide %in% c("V13Y10-020","V13B23-331"),]
dim(demo)
stopifnot(dim(demo)[1]==120)

## Build basic SPE
cat("\nBuild spe start time:", format(Sys.time(), tz="EST"), "\n")
round(system.time(spe <- read10xVisiumWrapper(
	demo$sample_path,
	demo$sample_id,
	type = "HDF5",
	data = "raw",
	images = c("lowres"),
	load = FALSE,
	reference_gtf = file.path("/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/","genes", "genes.gtf")
))/60,2)

#change colnames to be completely unique
spotcodes = paste(rownames(colData(spe)), spe$sample_id, sep="_")
colnames(spe) <- spotcodes
rownames(colData(spe)) <- spotcodes

## Complete colData
#merge info based on only shared colname (sample_id)
new_col <- merge(colData(spe), demo)
#remove unwanted columns
new_col <- new_col[,-grep("^10x", colnames(new_col))]
new_col <- new_col[,!colnames(new_col) %in% c("ManualAnnotation","sample_path")]
#make sure order is the same
new_col <- new_col[match(spe$key, new_col$key),]
stopifnot(identical(new_col$key, spe$key))
rownames(new_col) <- rownames(colData(spe))
colData(spe) <- new_col

cat("\n\nnew colData\n")
head(colData(spe))

cat("\n\nFinal spe\n")
for(i in reducedDimNames(spe)) 
  reducedDim(spe, i) <- NULL
spe

#save
#DelayedArray:::set_verbose_block_processing(TRUE)
#start.time1 = Sys.time()
#cat("\nStart HDF5Array save:"); start.time1
#cmtx = as(counts(spe), "HDF5Matrix")
#writeHDF5Array(cmtx, filepath="processed-data/02_build_spe/spe_n120_raw-counts_backup",
#               name="counts")
#cat("\nTime elapsed (writeHDF5Array):",
#        round(difftime(Sys.time(), start.time1, units="mins"),2), "minutes\n")


start.time2 = Sys.time()
cat("\nStart HDF5SummExp save:"); start.time2           
saveHDF5SummarizedExperiment(spe, dir="processed-data/02_build_spe", prefix="spe_n120_",
        #chunkdim=getHDF5DumpChunkDim(c(1,ncol(spe))),
        chunkdim=c(100,500),
	verbose=F)
cat("\nTime elapsed (saveHDF5SummarizedExperiment):",
        round(difftime(Sys.time(), start.time2, units="mins"),2), "minutes\n")

if(file.exists("spe_tracker_current.txt")) {
  if(file.exists("spe_tracker_archive.txt")) {
    x = readLines("spe_tracker_current.txt")
    write(c("########", paste("ARCHIVED",format(Sys.time(), tz="EST"),"EST"), "########"), "spe_tracker_archive.txt", append=TRUE)
  } else {
    writeLines(c("########", paste("ARCHIVED",format(Sys.time(), tz="EST"),"EST"), "########"), file("spe_tracker_archive.txt"))
    x = readLines("spe_tracker_current.txt")
    write(x, "spe_tracker_archive.txt", append=TRUE)
    close(file("spe_tracker_archive.txt"))
  }
}
writeLines(c(paste("Created raw spe (HDF5) on",format(Sys.time(), tz="EST"),"EST"),
             "New file location: processed-data/02_build_spe/spe_n120_",
             "Source code: code/02_build_spe/02_raw_spe_HDF5.r",
             "*","*","*"),
           file("spe_tracker_current.txt"))
close(file("spe_tracker_current.txt"))

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
