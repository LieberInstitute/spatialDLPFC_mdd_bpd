args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
    stop("Usage: Rscript 02_spe-from-GEO.R <geo_bundle> <sample_ids> <published_geo_dir>")
}
nest_dir <- args[[1L]]
sample_set <- read.table(args[[2L]], header=FALSE, stringsAsFactors=FALSE)[[1L]]
published_geo_dir <- args[[3L]]

suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
})

# BUILD BASIC SPE #
donor.data <- read.csv(file.path(nest_dir, "GEO_donor_data.csv"), row.names=1)

cat("\nBuild spe start time:", format(Sys.time(), tz="EST"), "\n")
round(system.time(spe <- read10xVisium(
	samples = file.path(nest_dir, sample_set, "outs"),
	sample_id = donor.data[sample_set, "sample_id"],
	type = "sparse",
	data = "raw",
	images = c("lowres"),
	load = FALSE
))/60,2)

#change colnames to be completely unique
spe$key = paste(rownames(colData(spe)), spe$sample_id, sep="_")
colnames(spe) <- spe$key
rownames(colData(spe)) <- spe$key

# then merge and re-order
cdata = merge(colData(spe), donor.data)
rownames(cdata) <- cdata$key

# to match colData columns and order with what is in obs supp file
#"sample_id"       "brnum"           "MBv_sample"     
#"diagnosis"       "sex"             "age"             "seq"            
#"key"             "array_row"       "array_col"       "in_tissue"
colData(spe) <- cdata[spe$key,c("sample_id","brnum","MBv_sample","dx","sex",
  "age","RIN",
  "batch","slide",
  "key","array_row","array_col","in_tissue")]

spe$sum_umi <- colSums(counts(spe))
spe$sum_gene <- colSums(counts(spe) > 0)
is_mito <- grep("^MT-", rowData(spe)$symbol)
spe$expr_chrM <- colSums(counts(spe)[is_mito, , drop = FALSE])
spe$expr_chrM_ratio <- spe$expr_chrM / spe$sum_umi


# LINK IMAGES TO THE PUBLISHED GEO BUNDLE #
## Without this section, the image path is to the cache/work dir which may be temporary or deleted later on
## This section replaces that path with the publish output dir path
image_data <- imgData(spe)
image_samples <- match(image_data$sample_id, donor.data$sample_id)
stopifnot(!anyNA(image_samples))
for (i in seq_len(nrow(image_data))) {
    image <- image_data$data[[i]]
    stopifnot(is(image, "StoredSpatialImage"))
    published_image <- file.path(published_geo_dir, rownames(donor.data)[image_samples[i]],
                                 "outs", "spatial", basename(imgSource(image)))
    # Publication is asynchronous: the final path need not exist inside this task.
    # Construct the stored reference without normalizing it against the work directory.
    image_data$data[[i]] <- new("StoredSpatialImage", path=published_image,
                                metadata=metadata(image))
}
imgData(spe) <- image_data

start.time2 = Sys.time()
cat("\nStart HDF5SummExp save:"); start.time2           
saveHDF5SummarizedExperiment(spe, ".",#dir="02_build_spe", 
    prefix=paste0("spe_n", length(sample_set), "_"),
    chunkdim=c(100,500),
	verbose=F)
cat("\nTime elapsed (saveHDF5SummarizedExperiment):",
        round(difftime(Sys.time(), start.time2, units="mins"),2), "minutes\n")
cat("\nFile prefix:", paste0("spe_n", length(sample_set), "_"), "\n")
