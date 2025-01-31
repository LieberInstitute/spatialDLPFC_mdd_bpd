setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(spdep)
	library(scater)
	library(dplyr)
	library(here)
})
set.seed(123)

start.time = Sys.time()
cat("Start time (loadHDF5SummarizedExperiment):", format(start.time),"\n")
spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_")
cat("Time elapsed (loadHDF5SummarizedExperiment):", round(difftime(Sys.time(), start.time, units="mins"),2), "minutes\n")

l1 = unique(spe$sample_id)
names(l1) = lapply(l1, function(x) unique(colData(spe)[spe$sample_id==x,"brain"]))
l1 = lapply(l1, function(x) colData(spe)[colData(spe)$sample_id==x,])

cat("\nLabel spots as off tissue, tissue perimeter, or image perimeter:", format(Sys.time()), "\n")
newColData = do.call(rbind, lapply(l1, function(x) {
	m1 = as.matrix(x[,c("array_row","array_col")])
	rownames(m1) = x$key
	nbors = dnearneigh(m1, d1=1, d2=2.1, bounds=c("GE","LE")) #d2 of 2.1 identifies closest 6 spots
	names(nbors) = x$key
	index= 1:nrow(x)
	off.tissue.index = index[x$in_tissue==FALSE]
	off.nbors = sapply(nbors, function(y) length(intersect(y, off.tissue.index)))
	x$keep = "keep"
	x[off.nbors>0,"keep"] = "tissue perimeter"
	x[off.tissue.index,"keep"] = "off tissue"

	east = max(x$array_row)
	west = min(x$array_row)
	north = max(x$array_col)
	south = min(x$array_col)
	cardinal.check = x$array_row %in% c(east, west) | x$array_col %in% c(north, south)
	x[x$keep!="off tissue" & cardinal.check,"keep"] = "image perimeter"
	return(x)
}))
Sys.time()
spe$keep_spots = newColData$keep

#now use 3MAD threshold to ID edges for removal, but first need to get rid of off tissue spots
spe2 = spe[,spe$keep_spots!="off tissue"]
spe2$lg10.umi = log10(spe2$sum_umi)
spe2$lg10.genes = log10(spe2$sum_gene)
cat("\n3MAD outlier to ID poor quality edges:", format(Sys.time()),"\n")
colData(spe2)$umi_3MAD.outlier_slide = isOutlier(spe2$lg10.umi, batch=spe2$slide, type="lower", nmads=3)
colData(spe2)$genes_3MAD.outlier_slide = isOutlier(spe2$lg10.genes, batch=spe2$slide, type="lower", nmads=3)

colData(spe2)$umi_3MAD.outlier_sample = isOutlier(spe2$lg10.umi, batch=spe2$sample_id, type="lower", nmads=3)
colData(spe2)$genes_3MAD.outlier_sample = isOutlier(spe2$lg10.genes, batch=spe2$sample_id, type="lower", nmads=3)

colData(spe2)$umi_3MAD.outlier = spe2$umi_3MAD.outlier_sample | spe2$umi_3MAD.outlier_slide
colData(spe2)$genes_3MAD.outlier = spe2$genes_3MAD.outlier_sample | spe2$genes_3MAD.outlier_slide
Sys.time()

#then ID samples for edge removal
l2 = unique(spe2$sample_id)
names(l2) = lapply(l2, function(x) unique(colData(spe2)[spe2$sample_id==x,"brain"]))
l2 = lapply(l2, function(x) colData(spe2)[colData(spe2)$sample_id==x,])

source(here("code","03_QC","01_edge_functions.r"))
cat("\nfindEdges() loop:", format(Sys.time()),"\n")
qualifying.edges = lapply(l2, function(x) 
	list("array_col"=findEdges(x, coord.dim="array_col"),
	"array_row"=findEdges(x, coord.dim="array_row")
	)
)
Sys.time()
outlierList = lapply(seq_along(qualifying.edges), function(x) {
	e.out = c()
	if(!is.na(qualifying.edges[[x]][[1]])) {
		e.out = c(e.out, idEdge(l2[[x]], qualifying.edges[[x]][1]))
	}
	if(!is.na(qualifying.edges[[x]][[2]])) {
		e.out = c(e.out, idEdge(l2[[x]], qualifying.edges[[x]][2]))
	}
	e.out
})

colData(spe)$edge_outlier = FALSE
colData(spe)[spe$key %in% unique(unlist(outlierList)),"edge_outlier"] = TRUE

table(colData(spe)[,c("keep_spots","edge_outlier")])

cat("Save colData:",format(Sys.time()),"\n")
write.csv(colData(spe), "processed-data/03_QC/spe_n120_edge-detection_colData.csv", row.names=T)

#save(spe, file=here("processed-data","03_QC","spe_demo.Rdata"))
#write(c(paste("*** Modified spe_demo on",format(Sys.time(), tz="UTC"),"UTC"),
#	paste("*** Old file location:",here("processed-data","02_build_spe","spe_demo.Rdata")),
#        paste("*** New file location:",here("processed-data","03_QC","spe_demo.Rdata")),
#        paste("*** Source code:",here("code","03_QC","01_revised_edge_detection.r")),
#        "*****","*****","*****"), here("spe_tracker_current.txt"), append=TRUE)

#cat("\n\nFiltering spe...\n")
#cat("Old dim:",dim(spe),"\n")
#spe = spe[,spe$keep_spots!="off tissue" & spe$edge_outlier==FALSE]
#cat("New dim:",dim(spe),"\n")
#cat("Saving filtered spe...\n")
#save(spe, file=here("processed-data","03_QC","spe_demo-filt.Rdata"))
#write(c(paste("***** Filtered spe_demo on",format(Sys.time(), tz="UTC"),"UTC"),
#        paste("***** Old file location:",here("processed-data","03_QC","spe_demo.Rdata")),
#        paste("***** New file location:",here("processed-data","03_QC","spe_demo-filt.Rdata")),
#        paste("***** Source code:",here("code","03_QC","01_revised_edge_detection.r")),
#        "*******","*******","*******"), here("spe_tracker_current.txt"), append=TRUE)

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
