args <- commandArgs(trailingOnly = TRUE)
sample_set <- args[[1L]]
maxseconds <- args[[2]]
raw_archive <- args[[3L]]
if (!file.exists(raw_archive)) {
  stop("GEO raw archive not found: ", raw_archive)
}
suppressPackageStartupMessages({
	library(GEOquery)
})

out_dir <- "GSE333974"
tmp_dir <- tempfile("geoq", tmpdir = tempdir())
create_dir <- function(path, recursive = FALSE) {
  if (!dir.create(path, recursive = recursive, showWarnings = FALSE)) {
    stop("Could not create directory (it must be new): ", path)
  }
}
create_dir(out_dir)
create_dir(tmp_dir)



# EXTRACT DONOR DATA FROM GEO MATRIX #

## Use direct URLs for the two platform-specific series matrices.
options(timeout = maxseconds) # increase timeout limit to 5 min for downloading large files
geo_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE333nnn/GSE333974"
download_geo_file <- function(subdir, filename) {
  url <- paste(geo_url, subdir, filename, sep = "/")
  destination <- file.path(tmp_dir, filename)
  if (download.file(url, destination, mode = "wb") != 0L) {
    stop("GEO download failed: ", url)
  }
  destination
}

series_matrices <- c("GSE333974-GPL24676_series_matrix.txt.gz",
                     "GSE333974-GPL34281_series_matrix.txt.gz")
metadata_columns <- c("supplementary_file_9", "instrument_model", "description",
                      "age:ch1", "agent:ch1", "batch:ch1", "Sex:ch1")
pd <- lapply(series_matrices, function(filename) {
  series_matrix <- download_geo_file("matrix", filename)
  metadata <- pData(getGEO(filename = series_matrix, getGPL = FALSE))
  metadata[,metadata_columns, drop=FALSE]
})
pd <- do.call(rbind, pd)

# Preserve GSM accession ordering before selecting the requested sample subset.
pd <- pd[order(rownames(pd)),,drop=FALSE]
colnames(pd) <- c("supplementary_file_9", "instrument_model", "RIN",
                  "age", "dx", "batch", "sex")
pd$RIN <- as.numeric(substr(pd$RIN, start=5, stop=7)) # RIN was rounded to 2 sigfig (0's were dropped)
pd$age <- as.numeric(pd$age)
pd$sex <- ifelse(pd$sex=="Male", "M", "F")


## Expand GSM name into component metadata.
sample_str <- sub("_matrix[.]mtx[.]gz$", "", basename(pd$supplementary_file_9))
# Keep sample_str first: the SPE build reads it as the donor table row names.
tmp = data.frame(sample_str=sample_str, pd, row.names=rownames(pd))

sdata = do.call(rbind.data.frame, strsplit(tmp$sample_str, "_"))

sdata = data.frame("GEO_sample"=sdata[,1], "MBv_sample"=paste(sdata[,2], sdata[,3], sep="_"),
  #"dx"=sdata[,4], "sex"=sdata[,5],
  "brnum"=sapply(strsplit(sdata[,6], "-"), function(x) x[[1]]),
  "slide"=substr(sdata[,6], start=8, stop=20), "array"=sdata[,7])
sdata$sample_id = paste(sdata$slide, sdata$array, sep="_")


## complete and save donor data
donor.data = cbind(tmp, sdata)
stopifnot(identical(rownames(donor.data), donor.data$GEO_sample))

write.csv(donor.data, file.path(out_dir, "GEO_donor_data.csv"), row.names=FALSE)




# EXTRACT SPACERANGER OUTPUTS FROM GEO #

## Extract the staged archive, then restructure the selected samples.
tar_dir <- file.path(tmp_dir, "extracted_raw")
stopifnot(untar(raw_archive, exdir = tar_dir) == 0L)

## then nest to match what read10xVisium expects spaceranger output to look like
matrix_dir <- "outs/raw_feature_bc_matrix"
spatial_dir <- "outs/spatial"
destinations <- c(
  "_barcodes.tsv.gz" = matrix_dir,
  "_features.tsv.gz" = matrix_dir,
  "_matrix.mtx.gz" = matrix_dir,
  "-scalefactors_json.json.gz" = spatial_dir,
  "-tissue_lowres_image.png.gz" = spatial_dir,
  "-tissue_positions.csv.gz" = spatial_dir
)

sids = donor.data$sample_str
nest_dir = out_dir
#nest_dir = file.path(out_dir, "nested")
#create_dir(nest_dir)

for (i in 1:sample_set) {
  sample_dir <- file.path(nest_dir, sids[i])
  create_dir(sample_dir)
  create_dir(file.path(sample_dir, matrix_dir), recursive = TRUE)
  create_dir(file.path(sample_dir, spatial_dir))
  
  f1 = paste0(file.path(tar_dir, sids[[i]]), names(destinations))
  stopifnot(sum(file.exists(f1))==6)
  
  # unzip spatials or else read10xVisium won't work
  for(j in f1[4:6]) R.utils::gunzip(j, destname = gsub("[.]gz$", "", j), remove = FALSE)
  
  # remove leading symbols for new file name
  target_path = file.path(sample_dir, destinations, substr(names(destinations), 2, 100))
  # remove .gz from spatial new file names
  target_path[4:6] = sapply(target_path[4:6], function(x) substr(x, 0, nchar(x)-3))

  copied <- file.copy(from= f1, to = target_path)
  stopifnot(all(copied))
}

# PRESERVE SAMPLE IDS #
write.table(sids[1:sample_set], file.path("GEO_sample_ids.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
