setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(Matrix)
})

# load full spe with two samples for Br5366 (n=120)
spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_")
dim(spe) # 36601 599034

# load spot-level QC file, copy QC inclusion
cdata = read.csv("processed-data/03_QC/colData_edges-problem-areas_spotsweeper_FINAL.csv", row.names=1)
stopifnot(identical(rownames(cdata), rownames(colData(spe))))
spe$remove_spots = cdata$remove_spots

# remove extra Br5366 sample with lower QC b/w the two
spe = spe[,spe$sample_id!="V13F27-338_C1"]
dim(spe) # 36601 594044

# tally number of spots with detected expression per gene
rowData(spe)$n_spots = rowSums(counts(spe)>0)

# extract the gene lists from processed pseudobulk objects
tmp <- readRDS("processed-data/16_globus/MBv_pseudobulk_domain-sp.rds")
rowData(spe)$SVG = rownames(spe) %in% rownames(tmp)[rowData(tmp)$SVG]
rowData(spe)$DE_domain.sp = rownames(spe) %in% rownames(tmp)[rowData(tmp)$DE_input]
tmp <- readRDS("processed-data/16_globus/MBv_pseudobulk_domain-ct.rds")
rowData(spe)$DE_domain.ct = rownames(spe) %in% rownames(tmp)[rowData(tmp)$DE_input]
tmp <- readRDS("processed-data/16_globus/MBv_pseudobulk_domain-ct-eQTL.rds")
rowData(spe)$grnboost2 = rownames(spe) %in% rownames(tmp)[rowData(tmp)$GRNBOOST2_input]

# extract colData columns
colData(spe) <- colData(spe)[,c("sample_id", "brnum", "MBv_sample", #ID keys
                                "condition", "sex", "age", #experimental group
				"seq", #sequencing batch
                                "key", "array_row", "array_col", "in_tissue", #spatial metrics
                                "sum_umi", "sum_gene", "expr_chrM_ratio", "remove_spots")] #QC metrics
# round age to tenths for privacy
colData(spe)$age = signif(colData(spe)$age, 3)
# update language
colnames(colData(spe))[4] = "diagnosis"
spe$diagnosis = gsub("BPD", "BD", spe$diagnosis)

# add spatial domain annotations
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
colData(spe)$domain_sp = NA
colData(spe)[rownames(cdata), "domain_sp"] = as.character(factor(cdata$smoothed_k9_1663_f, levels=c("L1","L2","L3/4","L5","L6","WM","low UMI","Vasc","GABA"),
                              labels=c("L1","L2","L3/4","L5","L6","WM","dropped.lowUMI","dropped.Vasc","dropped.GABA")))

res = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc30_red-precast-kweight-50-low-res.csv", row.names=1)
colData(spe)$domain_ct = NA
colData(spe)[rownames(res), "domain_ct"] = as.character(factor(res$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"),
                            labels=c("M/V","Ast","L2/3","L2/3","L4","Inb","L5","L6","Olg")))
#                domain_ct
#domain_sp          Ast   Inb  L2/3    L4    L5    L6   M/V   Olg
#  dropped.GABA      18   307     2     2     2     5     5     0
#  dropped.lowUMI  1225  3077   348  2961   650  2451  3254  8082
#  dropped.Vasc     140    57     1    22    40    95  4093    39
#  L1             33693  9609   163   333    59  6996 12695   851
#  L2             34497 16959 43035  3006   170  1512   912    80
#  L3/4           22824 15105 54999 31704 12928  5347   781   130
#  L5              3452  2439   105  1080 64091  1246   315   193
#  L6             12046  4899  1091  1812  2279 62378  1512  2037
#  WM              7304   288     7    63    84   658   952 29653

# add n109 indicator for SCENIC
npas4.outliers = c("Br5666","Br5594","Br5599","Br5448","Br6316","Br6021","Br5993","Br6192","Br5454","Br8073")
colData(spe)$SCENIC_n109 = !spe$brnum %in% npas4.outliers

# update spatial coords spot barcodes
rownames(spatialCoords(spe)) <- colnames(spe)

# write files
#write.csv(spatialCoords(spe), "processed-data/16_globus/MBv_n119_spatial-coords.csv", row.names=T)
##write.csv(rowData(spe), "processed-data/16_globus/MBv_n119_features.csv", row.names=T)
#write.csv(colData(spe), "processed-data/16_globus/MBv_n119_observations.csv", row.names=T)

#Sys.time()
#dgc_matrix <- as(counts(spe), "dgCMatrix")
#writeMM(dgc_matrix, file = "processed-data/16_globus/MBv_n119_raw-counts.mtx")
#cat("MM saved")

# subset to qc filtered object then re-calculate nspots detected per gene since this is what we used to find the 28,965 genes moving forward
spe = spe[,spe$remove_spots==FALSE]
dim(spe)
rowData(spe)$n_spots = rowSums(counts(spe)>0)
write.csv(rowData(spe), "processed-data/16_globus/MBv_n119_features.csv", row.names=T)

## Reproducibility information
cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
