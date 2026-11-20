args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
    stop("Usage: Rscript 05_actualize-QC-filters.R <assay_h5> <qc_coldata> <keep_generated_qc>")
}
h5Path <- args[[1L]]
fname = basename(h5Path)
h5_prefix = sub("assays\\.h5$", "", fname)

keep_generated_qc <- tolower(args[[3L]])
if (!keep_generated_qc %in% c("true", "false")) {
    stop("keep_generated_qc must be true or false.")
}
keep_generated_qc <- keep_generated_qc == "true"


# functions and setup for pulling in rowData from GEO
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


suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(HDF5Array)
})


# LOAD QC ANNOTATIONS AND ORIGINAL GSE333974 FILTERS #

qc.df <- read.csv(args[[2L]], row.names=1)
obs_path <- download_geo_file("suppl", "GSE333974_MBv_n119_observations.csv.gz")
pub.qc = read.csv(obs_path, row.names=1)
stopifnot(all(c("key", "sample_id", "remove_spots") %in% colnames(pub.qc)))

## Make CSV comparison agreeable regardless of subset or any ordering differences
#### Separate check for subset
pub.qc <- pub.qc[pub.qc$sample_id %in% unique(qc.df$sample_id), ]
#### Vs separate check for order
pub.qc <- pub.qc[match(qc.df$key, pub.qc$key), ]


# COMPARE GENERATED AND PUBLISHED FILTERS #

generated_remove_spots <- qc.df$remove_spots

comparison <- table(
    generated_remove_spots = factor(generated_remove_spots, levels=c(FALSE, TRUE)),
    published_remove_spots = factor(pub.qc$remove_spots, levels=c(FALSE, TRUE))
)
cat("\nGenerated versus published remove_spots:\n")
print(comparison)
write.csv(as.data.frame(comparison), "QC_filter_comparison.csv", row.names=FALSE)

different <- generated_remove_spots != pub.qc$remove_spots
slide <- sub("_[A-D][1-4]$", "", qc.df$sample_id)
differences <- data.frame(
    key=qc.df$key[different], sample_id=qc.df$sample_id[different],
    slide=slide[different], generated_remove_spots=generated_remove_spots[different],
    published_remove_spots=pub.qc$remove_spots[different]
)
write.csv(differences, "QC_filter_differences.csv", row.names=FALSE)
message(sum(different), " spots differ between generated and published QC filters.")
if (any(slide[different] != "V13F27-338")) {
    stop("QC differences occur outside the expected slide V13F27-338: ",
         paste(unique(slide[different & slide != "V13F27-338"]), collapse=", "))
}
if (any(different)) {
    message("All QC differences are restricted to V13F27-338. Published QC used four ",
            "samples on this slide; the regenerated QC uses three, changing the slide-level ",
            "3MAD thresholds and potentially downstream local outlier flags.")
}

qc_source <- if (keep_generated_qc) "generated" else "published"
if (keep_generated_qc) {
    warning("Keeping generated QC filters. This departs from the published dataset's QC filters.",
            call.=FALSE)
} else {
    message("Keeping published QC filters for identical reproduction of the published filtering.")
}
qc.df$generated_remove_spots <- generated_remove_spots
qc.df$published_remove_spots <- pub.qc$remove_spots
qc.df$remove_spots <- if (keep_generated_qc) generated_remove_spots else pub.qc$remove_spots
qc.df$qc_filter_source <- qc_source

write.csv(qc.df, "colData_QC.csv", row.names=TRUE)


# LOAD SPE, ADD QC ANNOTATIONS, AND SAVE #

spe <- loadHDF5SummarizedExperiment(dir=".", prefix=h5_prefix)
stopifnot(identical(colnames(spe), rownames(qc.df)))
for (column in colnames(qc.df)) {
    colData(spe)[[column]] <- qc.df[[column]]
}
metadata(spe)$QC_filter_source <- qc_source
metadata(spe)$QC_reference <- obs_path

quickResaveHDF5SummarizedExperiment(spe)

