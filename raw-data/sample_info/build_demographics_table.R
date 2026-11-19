setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')

# Compile demo table of all info used in any table or figure/SF panel in manuscript #

#### This contains too many potentially identifying characteristics to be able to be
#### made public as is (linking each metric to each donor). This information is 
#### presented at the aggregate / abstracted level when included in the manuscript

source("raw-data/sample_info/getMBvSampleInfo_function.R")

demo = getMBvSampleInfo(REDCapFile="Visium_DATA_2025-01-22_1406.csv",
    demoFile="DLPFC_cross-disorders_demographics_MBv.csv")

# remove sample missing spacerange (run with NAc i think)
demo = demo[demo$sample_id!="V13B23-339_A1",]

# remove re-run slides
demo = demo[!demo$slide %in% c("V13Y10-020","V13B23-331"),]

# remove repeated NTC M 
demo = demo[!demo$sample_id=="V13F27-338_C1",]

# clean up
demo$dx = gsub("BPD","BD", demo$condition)
demo$condition <- NULL
# 3 sig digits for privacy
demo$age = signif(demo$age, 3)
# seq/rnd expansion
demo$instrum = factor(demo$seq, levels = c("SC-TC", "SKCCC", "Psomagen-1", "Psomagen-2"),
    labels = c("Novaseq 6000", "Novaseq 6000", "Novaseq X", "Novaseq X"))
demo$batch = demo$seq
demo$seq <- NULL

# add additional variables
cdata = read.csv("raw-data/sample_info/Page_BMI_Smoking_CrossDisorder_120825.csv")
demo = merge(demo, cdata[,c("BrNum", "BMI", "Smoking")], by.x="brnum", by.y="BrNum")

df1 = read.csv("raw-data/sample_info/MBv_suicide.csv")
#demo = merge(demo, df1[,c("BrNum","Manner.Of.Death","Violent.method.")], 
#    by.x="brnum", by.y="BrNum", all.x=T)
demo$suicide = demo$brnum %in% df1$BrNum

df2 = read.csv("raw-data/sample_info/Bipolar_subtype.csv")
demo = merge(demo, df2[,c("BrNum","BP.subtype")], 
    by.x="brnum", by.y="BrNum", all.x=T)

# add geo title as row.names
rownames(demo) = paste(paste(demo$MBv_sample, demo$dx, demo$sex, demo$brnum, sep="_"), demo$sample_id, sep="-")
dim(demo) # 119 17

demo <- demo[sort(rownames(demo)),c(
    "sample_id", "brnum", "MBv_sample", #ID keys
    "dx", "sex", "BP.subtype", #experimental group
    "age","RIN","PMI","BMI","Smoking","suicide", #demo/death conditions
    "instrum","batch","round","slide","array")] #sequencing batch

write.csv(demo, "raw-data/sample_info/n119_demographics.csv", row.names=T)