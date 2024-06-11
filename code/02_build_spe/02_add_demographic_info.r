setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(dplyr)
	library(here)
})

d1 = read.csv(here("raw-data","sample_info","Big_240_DLPFC_Dissections.csv"))
d2 = read.csv(here("raw-data","sample_info","MDD_BPD_VisiumHE.csv"))
demo = left_join(d2, d1[,c("brain","age","sex","condition","PMI","RIN")], by="brain")

load(here::here("code", "REDCap", "REDCap_MBv.rda")) #don't know where this came from
demo = left_join(demo, REDCap_MBv[,c("slide","sample","array")], by=c("brain"="sample", "array"))

write.csv(demo, here("raw-data","sample_info","MBv_demographics_combined.csv"), row.names=FALSE)

load(here("processed-data","02_build_spe","spe_raw.Rdata"))

check.sample_id = setdiff(unique(spe$sample_id),unique(paste(demo$slide, demo$array, sep="_")))
if(length(check.sample_id)>0) {
cat("The following sample IDs were present in the spe object but not the demographic info:\n")
cat(check.sample_id)
cat("\n")
stop("Sample ID mismatch between demographics and spe object.")
}

mdata = mutate(as.data.frame(colData(spe)), slide=as.character(slide), array=as.character(array), brain=as.character(brnum))
mdata = left_join(mdata, demo, by=c("brain","slide","array"))
colnames(mdata)[grep("^array$",colnames(mdata))] = "position"

colData(spe) = cbind(colData(spe)[,c("key","sample_id")],
	mdata[,c("slide","position","brain","mbv_sample","age","sex","condition","PMI","RIN")],
	colData(spe)[,c("in_tissue","array_row","array_col",
	"sum_umi","sum_gene","expr_chrM","expr_chrM_ratio")])

save(spe, file=here("processed-data","02_build_spe","spe_demo.Rdata"))

write(c(paste("* Created spe_demo on",format(Sys.time(), tz="UTC"),"UTC"),
        paste("* Old file location:",here("processed-data","02_build_spe","spe_raw.Rdata")),
	paste("* New file location:",here("processed-data","02_build_spe","spe_demo.Rdata")),
        paste("* Source code:",here("code","02_build_spe","02_add_demographic_info.r")),
        "***","***","***"), here("spe_tracker_current.txt"), append=TRUE)

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
