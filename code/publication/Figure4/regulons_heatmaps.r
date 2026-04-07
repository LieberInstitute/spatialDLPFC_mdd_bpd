
library(dplyr)

aucell = read.csv("processed-data/10_SCENIC/spe-n109_13162-no-lowUMI_regulons-top20_AUCell.csv", row.names=1)
colnames(aucell)[1:(grep("seurat_label", colnames(aucell))-1)] = sapply(colnames(aucell)[1:(grep("seurat_label", colnames(aucell))-1)], 
	function(x) substr(x, start=0, stop=nchar(x)-3))

seurat_levels= c("M.V","Ast","L2.3","L4","Inb","L5","L6","Olg")
cond_sex = c("F NTC","F MDD", "F BPD","M NTC","M MDD","M BPD")

aucell$sample_id = substr(rownames(aucell), start=20, stop=50)
aucell$sample_id2 = paste(aucell$sample_id, aucell$seurat_label)
aucell$seurat_label= factor(aucell$seurat_label, levels=c("Micro.Vasc", "Astro", "L2.3", "L4", "Inhb", "L5", "L6", "Oligo"),
	labels=seurat_levels)
aucell$sex.group = paste(aucell$sex, aucell$condition)
aucell$clus_groups= factor(paste(aucell$sex.group, aucell$seurat_label),
	levels=as.character(outer(cond_sex, seurat_levels, paste)))

regulons <- read.csv("processed-data/10_SCENIC/spe-n109_13162-no-lowUMI_logcounts_regulons-top20_merged.csv", row.names=1)
regulons = regulons[regulons$set_size>5,]

# calculate averages per dx*sex*cluster
avg.df = group_by(aucell, clus_groups, sample_id) %>% summarise_at(regulons$TF, mean)

source("code/02_build_spe/getMBvSampleInfo_function.r")
demo = getMBvSampleInfo(REDCapFile="Visium_DATA_2025-01-22_1406.csv",
	demoFile="DLPFC_cross-disorders_demographics_MBv.csv")

cdata = read.csv("raw-data/sample_info/Page_BMI_Smoking_CrossDisorder_120825.csv")
demo = left_join(demo, cdata, by=c("brnum"="BrNum"))
avg.df = left_join(avg.df, demo[,c("sample_id","brnum","MBv_sample","age","RIN","BMI")], by="sample_id")

cor.mtx = cor(avg.df[,c(regulons$TF, "age", "RIN","BMI")])
cor.mtx[,c("age","RIN","BMI")]



#pc3
results_set = "seurat-pc30"
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-no-lowUMI_sample-seurat-pc30_norm-filt.Rdata")
colnames(colData(spe_pseudo))[grep("seurat", colnames(colData(spe_pseudo)))] = "cluster"

colnames(colData(spe_pseudo))[grep("subsets_mito_percent", colnames(colData(spe_pseudo)))] = "chrM_ratio"
spe_pseudo$pc3 = reducedDim(spe_pseudo, "PCA_1663")[,3]
spe_pseudo$sample_id2 = paste(spe_pseudo$sample_id, spe_pseudo$cluster)

cdata = read.csv("raw-data/sample_info/Page_BMI_Smoking_CrossDisorder_120825.csv")
#decide to use Smoking status
sdata = distinct(as.data.frame(colData(spe_pseudo)[,c("brnum","sample_id","sample_id2","condition","pc3","sex","age","RIN")]))
#dim(sdata)

sdata = left_join(sdata, cdata[,c("BrNum","BMI","Smoking")], by=c("brnum"="BrNum"))

avg.df = group_by(aucell, sample_id, seurat_label, sample_id2) %>% summarise_at(regulons$TF, mean) %>%
	left_join(sdata, by=c("sample_id","sample_id2"))
# 8 rows are NA b/c 10 pseudobulk samples were dropped
avg.df = avg.df[!is.na(avg.df$brnum),]
 
cor.mtx = cor(avg.df[,c(regulons$TF, "age","RIN","BMI","pc3")])

pheatmap(cor.mtx[regulons$TF,c("age","RIN","BMI","pc3")],
	breaks= seq(from = -1, to = 1, length.out = 101))

varList = levels(avg.df$seurat_label)
names(varList) = varList
corrList = lapply(varList, function(x) cor(avg.df[avg.df$seurat_label==x, c(regulons$TF, "age","RIN","BMI","pc3")]))

corrList = lapply(varList, function(x) cor(avg.df[avg.df$seurat_label==x, c(colnames(aucell)[1:15], "age","RIN","BMI","pc3")]))

hmpList = lapply(names(corrList), function(x) {
  tmp = pheatmap::pheatmap(corrList[[x]][colnames(aucell)[1:15],c("age","RIN","BMI","pc3")],
#  tmp = pheatmap::pheatmap(corrList[[x]][regulons$TF,c("age","RIN","BMI","pc3")], 
	 breaks= seq(from = -1, to = 1, length.out = 101),
	angle_col = 0, main=x, silent=T, fontsize=6,
	treeheight_row = 15, treeheight_col = 15)
  return(tmp[[4]])
})

do.call(grid.arrange, c(grobs=hmpList, ncol=4))


avg.mtx = as.matrix(avg.df[,-1])
rownames(avg.mtx) = avg.df$clus_groups
