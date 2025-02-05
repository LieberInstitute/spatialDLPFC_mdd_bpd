setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(dplyr)
	library(spdep)
	library(scater)
	library(SpotSweeper)
})
set.seed(123)

#load spe
spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_")
cat("Dim spe:",dim(spe),"\n")

#load coldata and select problem areas
cdata = read.csv("processed-data/03_QC/colData_edges-problem-areas.csv", row.names=1)
#all edges are appropriately identified except for from the following samples
### V13B23-302_A1 and V13B23-302_C1
cdata$true_edges = ifelse(cdata$slide=="V13B23-302", FALSE, cdata$edge_outlier_genes)

tmp = filter(cdata, in_tissue==TRUE, true_edges==FALSE) %>% mutate(lowumi = sum_umi<=100) %>%
  group_by(problem_areas_genes.id) %>%
  summarise(n_lowumi=sum(lowumi), n_spots=n(), prop_lowumi=n_lowumi/n_spots) %>%
  filter(n_spots>5, !is.na(problem_areas_genes.id))
remove.areas = unique(filter(tmp, prop_lowumi>=.5)$problem_areas_genes.id)
cat("Problem areas to remove based on >=50% spots with UMI<=100:", length(unique(remove.areas)),"\n")


#add edge info and subset
stopifnot(identical(rownames(cdata),rownames(colData(spe))))

spe$edge_outlier_genes = cdata$edge_outlier_genes
spe$true_edges = cdata$true_edges

spe$problem_areas_genes.id = cdata$problem_areas_genes.id
spe$problem_areas_genes.size = cdata$problem_areas_genes.size
spe$problem_areas_binary = spe$problem_areas_genes.id %in% remove.areas

spe$lowumi = spe$sum_umi<=100
spe$remove_spots = spe$problem_areas_binary | spe$true_edges | spe$lowumi

#copy removal vars to cdata for saving after spotsweeper
cdata$problem_aras_binary = spe$problem_areas_binary
cdata$lowumi = spe$lowumi
cdata$remove_spots = spe$remove_spots
spe = spe[,spe$in_tissue]
cat("Dim spe (in tissue):",dim(spe),"\n")

cat("Criteria for spot removal\n")
table(colData(spe)[,c("problem_areas_binary","lowumi","true_edges")])

cat("Total spots for removal:")
table(spe$remove_spots)

spe = spe[,spe$remove_spots==FALSE]
cat("Dim spe (after filtering):",dim(spe),"\n")

#run spotsweeper
cat(format(Sys.time()), "Calculate local outliers (umi counts)...","\n")
spe <- localOutliers(spe, metric = "sum_umi", direction = "lower", log = TRUE)
#colData(spe)[,c(ncol(colData(spe)),ncol(colData(spe))-2)] = NULL
colnames(colData(spe))[ncol(colData(spe))-1] = "umi_local.outlier"

cat(format(Sys.time()), "Calculate local outliers (n genes)...","\n")
spe <- localOutliers(spe, metric = "sum_gene", direction = "lower", log = TRUE)
#colData(spe)[,c(ncol(colData(spe)),ncol(colData(spe))-2)] = NULL
colnames(colData(spe))[ncol(colData(spe))-1] = "genes_local.outlier"

cat(format(Sys.time()), "Calculate local outliers (mito %)...","\n")
spe <- localOutliers(spe, metric = "expr_chrM_ratio", direction = "higher", log = FALSE)
#colData(spe)[,c(ncol(colData(spe)),ncol(colData(spe))-2)] = NULL
colnames(colData(spe))[ncol(colData(spe))-1] = "chrM.ratio_local.outlier"

#update coldata
true.spots = cdata$remove_spots==FALSE & cdata$in_tissue==TRUE
stopifnot(identical(rownames(cdata)[true.spots], rownames(colData(spe))))
cdata$umi_local.outlier = NA
cdata[true.spots, "umi_local.outlier"] = colData(spe)$umi_local.outlier
cdata$genes_local.outlier = NA
cdata[true.spots, "genes_local.outlier"] = colData(spe)$genes_local.outlier
cdata$chrM.ratio_local.outlier = NA
cdata[true.spots, "chrM.ratio_local.outlier"] =	colData(spe)$chrM.ratio_local.outlier

#save coldata
write.csv(cdata, "processed-data/03_QC/colData_edges-problem-areas_spotsweeper.csv", row.names=T)
cat("\nUpdated colData saved to: processed-data/03_QC/colData_edges-problem-areas_spotsweeper.csv\n")
#save(spe, file=here("processed-data","03_QC","spe_demo-filt.Rdata"))

#write(c(paste("******* Modified spe_demo-filt on",format(Sys.time(), tz="UTC"),"UTC"),
#        paste("******* Old location:",here("processed-data","03_QC","spe_demo-filt.Rdata")),
#	paste("******* New location:",here("processed-data","03_QC","spe_demo-filt.Rdata")),
#        paste("******* Source code:",here("code","03_QC","02_outliers.r")),
#        "*********","*********","*********"), here("spe_tracker_current.txt"), append=TRUE)
#
## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
