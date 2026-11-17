args <- commandArgs(trailingOnly = TRUE)
maxseconds <- args[[1L]]
h5Path <- args[[2L]]

fname = basename(h5Path)
h5_prefix = gsub("assays.h5$", "", fname)

suppressPackageStartupMessages({
	library("GEOquery")
	library("SpatialExperiment")
	library("HDF5Array")
    #library("rtracklayer") # will need when this is finished but for now just using the provided GSE333974 feature data
})

spe <- loadHDF5SummarizedExperiment(dir = ".", prefix = h5_prefix)

# load in rowData
## using reference GTF file
#gtf <- rtracklayer::import(reference_gtf)
## this is only the first line of many: https://github.com/LieberInstitute/spatialLIBD/blob/devel/R/read10xVisiumWrapper.R#L105

# if pulling from provided GSE333974 feature data
options(timeout = maxseconds) # increase timeout limit to 5 min for downloading large files
tmp_dir <- tempfile("geoq", tmpdir = tempdir())
dir.create(tmp_dir, showWarnings = FALSE)
geo_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE333nnn/GSE333974"
download_geo_file <- function(subdir, filename) {
  url <- paste(geo_url, subdir, filename, sep = "/")
  destination <- file.path(tmp_dir, filename)
  if (download.file(url, destination, mode = "wb") != 0L) {
    stop("GEO download failed: ", url)
  }
  destination
}
features_path <- download_geo_file("suppl", "GSE333974_MBv_n119_features.csv.gz")

fdata = read.csv(features_path, row.names=1)
stopifnot(identical(rownames(spe), rownames(fdata)))
rowData(spe) = as(fdata[,1:7], "DataFrame")


rowData(spe)

quickResaveHDF5SummarizedExperiment(spe)