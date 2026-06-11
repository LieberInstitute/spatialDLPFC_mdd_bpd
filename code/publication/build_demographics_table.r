setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(dplyr)
})

source("code/02_build_spe/getMBvSampleInfo_function.r")

demo = getMBvSampleInfo(REDCapFile="Visium_DATA_2025-01-22_1406.csv",
                        demoFile="DLPFC_cross-disorders_demographics_MBv.csv")

# remove sample missing spacerange (run with NAc i think)
demo = demo[demo$sample_id!="V13B23-339_A1",]

# remove re-run slides
demo = demo[!demo$slide %in% c("V13Y10-020","V13B23-331"),]

# remove repeated NTC M 
demo = demo[!demo$sample_id=="V13F27-338_C1",]

# add additional variables
cdata = read.csv("raw-data/sample_info/Page_BMI_Smoking_CrossDisorder_120825.csv")
cdata = select(cdata, brnum=BrNum, BMI, Smoking)

# merge and add title
demo = merge(demo, cdata)
demo$condition = gsub("BPD","BD", demo$condition)
demo$title = paste(demo$MBv_sample, demo$condition, demo$sex, demo$brnum, sep="_")

# 3 sig digits for privacy
demo$age = signif(demo$age, 3)

# order and save
demo = arrange(demo[,c("title","sample_id","brnum","MBv_sample","seq",
	"condition","sex","age","RIN","BMI","Smoking")], MBv_sample)

head(demo)

write.csv(demo, "processed-data/publication/demographics.csv", row.names=F)

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
