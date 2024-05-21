library(SpatialExperiment)
library(dplyr)
library(here)

load(here("processed-data","02_build_spe","spe_raw.Rdata"))
mdata = cbind.data.frame(colData(spe),
do.call(rbind, lapply(strsplit(spe$sample_id, split="_"), function(x) cbind.data.frame("slide"=x[[1]],"position"=x[[2]]))))

d1 = read.csv(here("raw-data","sample_info","all_brain_demographics_2024-05-21.csv"))
d2 = read.csv(here("raw-data","sample_info","mdd_bpd_brains_2024-05-21.csv"))
demo = left_join(d2, d1[,c("brain","sex","condition")], by=c("brain"))
mdata = left_join(mdata, demo, by=c("slide","position"))
 
colData(spe) = cbind(colData(spe)[,c("key","sample_id")],
	mdata[,c("slide","position","brain","sex","condition")],
	colData(spe)[,c("in_tissue","array_row","array_col",
	"sum_umi","sum_gene","expr_chrM","expr_chrM_ratio")])

save(spe, file=here("processed-data","02_build_spe","spe_demo.Rdata"))
