library(SpatialExperiment)
library(ggspavis)
library(here)

load(here("processed-data","03_QC","spe_demo-filt.Rdata"))

table(colData(spe)[,c("brain","exclude_east_edge")])

spe$exclude_outliers = spe$sum_local.outlier | spe$genes_local.outlier | spe$chrM.ratio_local.outlier | spe$exclude_east_edge
spe$extra_3MAD_outliers = spe$sum_3MAD.outlier_sample | spe$genes_3MAD.outlier_sample
spe$plot_outliers = factor(paste(spe$exclude_outliers, spe$extra_3MAD_outliers),
		levels=c("FALSE FALSE","FALSE TRUE","TRUE FALSE","TRUE TRUE"),
		labels=c("keep","3MAD flag","exclude","exclude"))

l1 = unique(spe$sample_id)
names(l1) = lapply(l1, function(x) unique(colData(spe)[spe$sample_id==x,"brain"]))
l1 = lapply(l1, function(x) colData(spe)$sample_id==x)

ggsave(plotSpots(spe, annotate = "plot_outliers", sample_id="sample_id",pal=c("lightgrey", "black", "red"), point_size=.2), 
filename=here("plots","03_QC","edge_detection_results.pdf"), width=12, height=12)
