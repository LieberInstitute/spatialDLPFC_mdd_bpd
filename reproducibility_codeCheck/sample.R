suppressPackageStartupMessages({
	library("GEOquery")
	library("SpatialExperiment")
	library("spatialLIBD")
	library("rtracklayer")
	library("HDF5Array")
})

#base_dir = "MBv_pipeline_revision"
getwd()

tmp_dir <- tempfile("geoquery-", tmpdir = tempdir())
dir.create(tmp_dir)

gse <- getGEO("GSE333974", destdir = tmp_dir)
#Found 1 file(s)
#GSE333974_series_matrix.txt.gz
pd = pData(gse[[1]])
tmp = pd[,44:46]
tmp[,1] = sapply(strsplit(sapply(strsplit(tmp[,1], "suppl\\/"), function(x) x[[2]]), "_matrix"), function(y) y[[1]])
colnames(tmp) = c("sample_str","dx","sex")

sdata = do.call(rbind.data.frame, strsplit(tmp$sample_str, "_"))

sdata = data.frame("GEO_sample"=sdata[,1], "MBv_sample"=paste(sdata[,2], sdata[,3], sep="_"),
           #"dx"=sdata[,4], "sex"=sdata[,5],
           "brnum"=sapply(strsplit(sdata[,6], "-"), function(x) x[[1]]),
           "slide"=substr(sdata[,6], start=8, stop=20), "array"=sdata[,7])
sdata$sample_id = paste(sdata$slide, sdata$array, sep="_")

donor.data = cbind(tmp, sdata)
stopifnot(identical(rownames(donor.data), donor.data$GEO_sample))

# can merge this with ST1 for age and rin demo data



supp_files <- getGEOSuppFiles('GSE333974', fetch_files = FALSE)
supp_files
# use first 4 files to build spe

# use last tar file to build from scratch spaceranger output
options(timeout = 300) # increase timeout limit to 5 min
supp_files <- getGEOSuppFiles('GSE333974', baseDir = tmp_dir, filter_regex = "RAW\\.tar")

# TO-DO: FIX PATHS
# figure out where GEOquery saved or predetermine that to be to tmp/cache is
# the extracted raw can also be in tmp path
tar_dir <- file.path(tmp_dir, "GSE333974/extracted_raw")
untar(file.path(tmp_dir, "GSE333974/GSE333974_RAW.tar"), exdir = tar_dir)
# untar creates flat dir

# nest files
create_dir <- function(path, recursive = FALSE) {
  if (!dir.create(path, recursive = recursive, showWarnings = FALSE)) {
    stop("Could not create directory (it must be new): ", path)
  }
}

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
nest_dir = file.path(tmp_dir, "GSE333974/nested")
create_dir(nest_dir)

for (i in 1:3) {
#for (i in seq_along(sids)) {
  sample_dir <- file.path(nest_dir, sids[i])
  create_dir(sample_dir)
  create_dir(file.path(sample_dir, matrix_dir), recursive = TRUE)
  create_dir(file.path(sample_dir, spatial_dir))
  
  f1 = paste0(file.path(tar_dir, sids[[i]]), names(destinations))
  stopifnot(sum(file.exists(f1))==6)
  
  # unzip spatials or else read10xVisium won't work
  for(j in f1[4:6]) gunzip(j, destname = gsub("[.]gz$", "", j), remove = FALSE)
  
  # remove leading symbols for new file name
  target_path = file.path(sample_dir, destinations, substr(names(destinations), 2, 100))
  # remove .gz from spatial new file names
  target_path[4:6] = sapply(target_path[4:6], function(x) substr(x, 0, nchar(x)-3))

  copied <- file.copy(from= f1, to = target_path)
  stopifnot(all(copied))
}

# build basic SPE
cat("\nBuild spe start time:", format(Sys.time(), tz="EST"), "\n")
round(system.time(spe <- read10xVisium(
	samples = file.path(nest_dir, sids[1:3], "outs"),
	sample_id = donor.data$sample_id[1:3],
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
                                #"age","seq",
                                "key","array_row","array_col","in_tissue")]


# load in rowData
## using reference GTF file
#gtf <- rtracklayer::import(reference_gtf)
## this is only the first line of many: https://github.com/LieberInstitute/spatialLIBD/blob/devel/R/read10xVisiumWrapper.R#L105

# if pulling from provided GSE333974 feature data
sf <- getGEOSuppFiles("GSE333974", baseDir = file.path(tmp_dir, "GSE333974"), makeDirectory = F,
	filter_regex = "feature")
fdata = read.csv(rownames(sf), row.names=1)
stopifnot(identical(rownames(spe), rownames(fdata)))
rowData(spe) = as(fdata[,1:7], "DataFrame")

spe$sum_umi <- colSums(counts(spe))
spe$sum_gene <- colSums(counts(spe) > 0)
is_mito <- grep("^MT-", rowData(spe)$gene_name)
spe$expr_chrM <- colSums(counts(spe)[is_mito, , drop = FALSE])
spe$expr_chrM_ratio <- spe$expr_chrM / spe$sum_umi


start.time2 = Sys.time()
cat("\nStart HDF5SummExp save:"); start.time2           
saveHDF5SummarizedExperiment(spe, dir="test_results/local_docker/", prefix="spe_n3_",
        #chunkdim=getHDF5DumpChunkDim(c(1,ncol(spe))),
        chunkdim=c(100,500),
	verbose=F)
cat("\nTime elapsed (saveHDF5SummarizedExperiment):",
        round(difftime(Sys.time(), start.time2, units="mins"),2), "minutes\n")
cat("\nFile path: test_results/local_docker/\nFile prefix: spe_n3_\n")


unlink(tmp_dir, recursive = TRUE)

cat("\n\n*** TEST END ***\n\n")
sessioninfo::session_info()
