#sensitivity analysis 
#covariates and colinearity

library(edgeR)
library(dplyr)
library(SpatialExperiment)
library(ggplot2)

cpList = readRDS("plots/colorPalettes.rds")

load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
spe_pseudo$pc3 = reducedDim(spe_pseudo, "PCA_1663")[,3]
comp_names = c("L1","L2","L3dot4","L5","L6","WM")
names(comp_names) = c("L1","L2","L3.4","L5","L6","WM")

cdata = read.csv("raw-data/sample_info/Page_BMI_Smoking_CrossDisorder_120825.csv")
colSums(is.na(cdata))

cdata = filter(cdata, BrNum %in% spe_pseudo$brnum)
table(cdata[,c("Nicotine","Cotinine")])
#.............Cotinine
#Nicotine     Negative Not Tested Positive
#..Negative         54          0        8 #8 negative for nicotine and positive for cotinine
#..Not Tested        0         37        0
#..Positive          0          0       20

table(cdata[,c("Smoking","Cotinine")])
#.............Cotinine
#Nicotine     Negative Not Tested Positive
#..Negative         54          0        8
#..Not Tested        0         37        0
#..Positive          0          0       20

#decide to use Smoking status
sdata = distinct(as.data.frame(colData(spe_pseudo)[,c("brnum","sample_id","condition","sex","age","PMI","RIN")]))
dim(sdata)

sdata = left_join(sdata, cdata[,c("BrNum","BMI","Smoking")], by=c("brnum"="BrNum"))
colSums(is.na(sdata))

#co-linearity? not an issue
cor(sdata[,c("age","PMI","RIN","BMI")])
#.....age         PMI         RIN         BMI
#age  1.00000000  0.17342499  0.05166935 -0.07249217
#PMI  0.17342499  1.00000000  0.06303647 -0.16666976
#RIN  0.05166935  0.06303647  1.00000000 -0.16601974
#BMI -0.07249217 -0.16666976 -0.16601974  1.00000000

table(sdata[,c("condition","sex","Smoking")])
#, , Smoking = No
#...........sex
#condition  F  M
#......NTC 17 17
#......MDD  7  7
#......BPD  4  8

#, , Smoking = Yes
#...........sex
#condition  F  M
#......NTC  3  3
#......MDD 13 12
#......BPD 16 12


#heatmap and percent variance
new.cdata = merge(colData(spe_pseudo), sdata, sort=F)
identical(spe_pseudo$total, new.cdata$total)
colData(spe_pseudo) <- new.cdata

exp.vars = c("smoothed_k9_1663","sample_id",
             "slide","seq",
             "condition","sex",
             "age","BMI","Smoking",
             "PMI","RIN",
             "sum","detected","nspots")
var.m = scater::getVarianceExplained(spe_pseudo, variables=c(exp.vars,"pc3"), 
                                     exprs_values="logcounts")
#decidedly not normal distribution
cor.var.m = cor(var.m, method="spearman")
col_annot = data.frame(colMeans(var.m))
colnames(col_annot) = "percVar"
round(col_annot, 2)
#...................percVar
#smoothed_k9_1663   18.56
#sample_id          28.55
#slide              13.16
#seq                 3.20
#condition           0.56
#sex                 0.34
#age                 0.55
#BMI                 0.29
#Smoking             0.28
#PMI                 0.38
#RIN                 0.30
#sum                 5.06
#detected            8.83
#nspots              4.41
#pc3                 6.58

ann_colors = list(
  percVar = colorRampPalette(c("white", "purple3", "black"), bias=1)(10)
)
hmp = pheatmap::pheatmap(cor.var.m,
                         annotation_col = col_annot, annotation_colors = ann_colors,
                         annotation_names_col=FALSE, annotation_legend=T)
plot(hmp[[4]])
