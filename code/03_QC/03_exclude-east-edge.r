setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(Seurat)
	library(here)
})

load(here("processed-data","03_QC","spe_demo-filt.Rdata"))

l1 = unique(spe$sample_id)
names(l1) = lapply(l1, function(x) unique(colData(spe)[spe$sample_id==x,"brain"]))
l1 = lapply(l1, function(x) spe[,colData(spe)$sample_id==x])

# step 1: find the eastern-most coordinate for all images
eastern.edge = lapply(l1, function(x) max(colData(x)[,"array_row"]))
n.eastern.edge = lapply(seq_along(l1), function(x) colData(l1[[x]])[l1[[x]]$array_row==eastern.edge[[x]],"array_row"])
names(n.eastern.edge) = names(eastern.edge)

# step 2: identify images that qualify for eastern edge exclusion
qualify = lapply(seq_along(l1), function(x) {
  subset1 = colData(l1[[x]])[l1[[x]]$array_row==eastern.edge[[x]],]
  perc = max(sum(subset1$sum_3MAD.outlier_sample)/nrow(subset1), sum(subset1$genes_3MAD.outlier_sample)/nrow(subset1))
  if(perc>.8 & length(n.eastern.edge[[x]])>20) {return(TRUE)}
  else {return(FALSE)}
})
names(qualify) = names(eastern.edge)

# step 3: generate adjacency matrix for qualifying images
l2 = l1[names(qualify)[unlist(qualify)]]

adj.mtx = lapply(l2, function(x) {
	count <- counts(x)
	a1 <- CreateAssayObject(count, assay = "RNA", min.features = 0, min.cells = 0)
	colData(x)$row = colData(x)$array_row
	colData(x)$col = colData(x)$array_col
	s1 = CreateSeuratObject(a1, meta.data = as.data.frame(colData(x)))
	a.mtx = DR.SC::getAdj(s1, platform="Visium")
	colnames(a.mtx) = colnames(count)
	rownames(a.mtx) = colnames(count)
	a.mtx
})

# step 4: step-wise exclude all outlier spots adjacent to eastern edge outliers
final.list = lapply(seq_along(l2), function(x) {
	cat(names(l2)[x],"\n")
	is.outlier = l2[[x]]$sum_3MAD.outlier_sample | l2[[x]]$genes_3MAD.outlier_sample
	out.mtx = adj.mtx[[x]][,is.outlier]
	new.seed = l2[[x]]$array_row==eastern.edge[[names(l2)[x]]] & is.outlier
	exclude.spots = colnames(adj.mtx[[x]])[new.seed]

	while(length(new.seed)>1) {
		r1 = colSums(out.mtx[new.seed,])>0
		new.seed = setdiff(colnames(out.mtx)[r1], exclude.spots)
		exclude.spots = union(new.seed, exclude.spots)
	}
	colData(l2[[x]])[exclude.spots,"key"] 
})

# step 5: create a new coldata column
spe$exclude_east_edge = FALSE
colData(spe)[spe$key %in% unlist(final.list),"exclude_east_edge"] = TRUE

save(spe, file=here("processed-data","03_QC","spe_demo-filt.Rdata"))

write(c(paste("Modified spe_demo-filt on",format(Sys.time(), tz="UTC"),"UTC"),
        paste("File location:",here("processed-data","03_QC","spe_demo-filt.Rdata")),
        paste("Source code:",here("code","03_QC","03_exclude-east-edge.r")),
        "*","*","*"), here("spe_tracker_current.txt"), append=TRUE)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
