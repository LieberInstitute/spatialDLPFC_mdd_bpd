args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
    stop("Usage: Rscript 01_find_edges-outlier-areas.R <assay_h5> <edge_functions>")
}
h5Path <- args[[1L]]
fname = basename(h5Path)
h5_prefix = sub("assays\\.h5$", "", fname)

# NEXTFLOW ENTRY REQUREMENT CHECK
required_files <- paste0(h5_prefix, c("assays.h5", "se.rds"))
if (!all(file.exists(required_files))) {
    stop("Missing SPE files: ", paste(required_files[!file.exists(required_files)], collapse = ", "))
}

suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(scuttle)
	library(raster)
	library(dplyr)
})

# LOAD SPE AND RASTER FUNCTIONS #

source(args[[2L]])
spe <- loadHDF5SummarizedExperiment(dir=".", prefix=h5_prefix)

## THIS WILL BE ABLE TO GO AFTER GEO SAMPLE MATRIX FILE IS UPDATED
# The GEO stage 1 SPE retains sample_id but omits the original slide column.
if (!"slide" %in% colnames(colData(spe)) && "sample_id" %in% colnames(colData(spe))) {
    if (!all(grepl("^V[^_]+_[A-D][1-4]$", spe$sample_id))) {
        stop("Cannot derive slide from sample_id. Supply an SPE with a slide column.")
    }
    spe$slide <- sub("_[A-D][1-4]$", "", spe$sample_id)
}
## END PART THAT CAN BE REMOVED

## Since starting a new nextflow entry point, check to make sure all of the expected column generated in Stage 1 are present
required_columns <- c("sum_umi", "sum_gene", "in_tissue", "slide", "sample_id", "array_row", "array_col")
if (!all(required_columns %in% colnames(colData(spe)))) {
    stop("SPE is missing required colData columns: ",
         paste(setdiff(required_columns, colnames(colData(spe))), collapse = ", "))
}




# 3MAD TO ID POOR QUALITY AREAS #

spe$lg10.genes = log10(spe$sum_gene)

## Based on MAD calculated per-slide
colData(spe)$genes_3MAD.outlier_slide = isOutlier(spe$lg10.genes, subset=spe$in_tissue, batch=spe$slide, type="lower", nmads=3)
## Based on MAD calculated per-sample
colData(spe)$genes_3MAD.outlier_sample = isOutlier(spe$lg10.genes, subset=spe$in_tissue, batch=spe$sample_id, type="lower", nmads=3)
# Union per-slide and per-sample outliers
colData(spe)$genes_3MAD.outlier_binary = colData(spe)$genes_3MAD.outlier_slide | colData(spe)$genes_3MAD.outlier_sample
colData(spe)$genes_3MAD.outlier_binary = ifelse(colData(spe)$in_tissue==FALSE, FALSE, colData(spe)$genes_3MAD.outlier_binary)


# PERFORM RASTER-BASED CONTIGUOUS AREA ID #

sampleList = unique(spe$sample_id)
names(sampleList) <- sampleList

## Edge detection
genes_edges = lapply(sampleList, function(x) {
	tmp = colData(spe)[spe$sample_id==x,c("in_tissue","array_row","array_col","genes_3MAD.outlier_binary")]
	clumpEdges(tmp[,-1], rownames(tmp)[tmp$in_tissue==FALSE])
})

cat("\nNumber of samples with edges detected")
table(sapply(genes_edges, length)>0) 
cat("\nSamples with edges detected\n")
names(genes_edges)[sapply(genes_edges, length)>0]

colData(spe)$edge_outlier = FALSE
colData(spe)[unlist(genes_edges),"edge_outlier"] = TRUE

## Contiguous areas of very low outliers
genes_probs = lapply(sampleList, function(x) {
	tmp = colData(spe)[spe$sample_id==x,c("in_tissue","array_row","array_col","genes_3MAD.outlier_binary")]
	problemAreas(tmp[,-1], 
		rownames(tmp)[tmp$in_tissue==FALSE], uniqueIdentifier=x, shifted=F)
})
genes_probs = do.call(rbind, genes_probs)
colData(spe)$MAD_area_id = NA
colData(spe)[genes_probs$spotcode,"MAD_area_id"] = genes_probs$clumpID
colData(spe)$MAD_area_size = 0
colData(spe)[genes_probs$spotcode,"MAD_area_size"] = genes_probs$clumpSize

## Restrict 3MAD outlier areas to those that are >20 spots and have >50% very low UMI
tmp = filter(as.data.frame(colData(spe)), in_tissue==TRUE) %>% mutate(lowumi = sum_umi<=100) %>%
  group_by(MAD_area_id) %>%
  summarise(n_lowumi=sum(lowumi), n_spots=n(), prop_lowumi=n_lowumi/n_spots) %>%
  filter(n_spots>20, !is.na(MAD_area_id))
remove.areas = unique(filter(tmp, prop_lowumi>=.5)$MAD_area_id)

colData(spe)$outlier_area = FALSE
colData(spe)[spe$MAD_area_id %in% remove.areas, "outlier_area"] = TRUE


# SAVE 3MAD-BASED OUTLIER DETECTED #
write.csv(colData(spe), "colData_QC.csv", row.names=T)

