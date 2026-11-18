args <- commandArgs(trailingOnly = TRUE)
maxseconds <- args[[1L]]
h5Path <- args[[2L]]
reference_gtf <- if (length(args) > 2L && nzchar(args[[3L]])) args[[3L]] else NA_character_

fname = basename(h5Path)
h5_prefix = gsub("assays.h5$", "", fname)

suppressPackageStartupMessages({
	library("GEOquery")
	library("SpatialExperiment")
	library("HDF5Array")
  library("rtracklayer")
})

spe <- loadHDF5SummarizedExperiment(dir = ".", prefix = h5_prefix)

# functions and setup for pulling in rowData from GEO
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

rowdata_from_GEO <- function(spe) {
  features_path <- download_geo_file("suppl", "GSE333974_MBv_n119_features.csv.gz")
  fdata = read.csv(features_path, row.names=1)
  stopifnot(identical(rownames(spe), rownames(fdata)))
  rowData(spe) = as(fdata[,1:7], "DataFrame")
  return(spe)
}



# UPDATE ROWDATA FROM GTF, ELSE FROM GEO #

if (!is.na(reference_gtf)) {
  
  # spatialLIBD code: https://github.com/LieberInstitute/spatialLIBD/blob/devel/R/read10xVisiumWrapper.R#L105

  if(file.exists(reference_gtf)) {
    gtf <- rtracklayer::import(reference_gtf)
    gtf <- gtf[gtf$type == "gene"]
    names(gtf) <- gtf$gene_id

    ## Match the genes
    message("Adding gene information from ", reference_gtf)
    match_genes <- match(rownames(spe), gtf$gene_id)

    if (all(is.na(match_genes))) {
    
        ## Protect against scenario where one set has GENCODE IDs and the other one has ENSEMBL IDs.
        warning("Gene IDs did not match. This typically happens when you are not using the same GTF file as the one that was used by SpaceRanger. For example, one file uses GENCODE IDs and the other one ENSEMBL IDs. read10xVisiumWrapper() will try to convert them to ENSEMBL IDs.", call. = FALSE)
        match_genes <- match(gsub("\\..*", "", rownames(spe)), gsub("\\..*", "", gtf$gene_id))
    
    }

    if (any(is.na(match_genes))) {
        
        warning("Dropping ", sum(is.na(match_genes)), " out of ", length(match_genes), " genes for which we don't have information on the reference GTF file. This typically happens when you are not using the same GTF file as the one that was used by SpaceRanger.", call. = FALSE)
        ## Drop the few genes for which we don't have information
        spe <- spe[!is.na(match_genes), ]
        match_genes <- match_genes[!is.na(match_genes)]
    
    }

    ## Keep only some columns from the gtf
    gtf_cols <- c("source", "type", "gene_id", "gene_version", "gene_name", "gene_type")
    mcols(gtf) <- mcols(gtf)[, gtf_cols[gtf_cols %in% colnames(mcols(gtf))]]

    ## Add the gene info to our SPE object
    rowRanges(spe) <- gtf[match_genes]

  } else {
    
    warning("Provided GTF file does not exist: ", reference_gtf,
            ". Using the provided GSE333974 feature data instead.", call. = FALSE)
    spe <- rowdata_from_GEO(spe)
  
  }
    
} else {

    warning("No reference GTF file provided. Using the provided GSE333974 feature data instead.", call. = FALSE)
    spe <- rowdata_from_GEO(spe)

}

quickResaveHDF5SummarizedExperiment(spe)
