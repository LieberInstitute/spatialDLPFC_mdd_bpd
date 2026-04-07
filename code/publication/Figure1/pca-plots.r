setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(scater)
  library(gridExtra)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(ggrastr)
})

set.seed(123)


cpList = readRDS("plots/colorPalettes.rds")

#PRECAST smoothed
load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")

p1 <- plotReducedDim(spe_pseudo, dimred="PCA_1663", ncomponents=2, colour_by = "smoothed_k9_1663", point_alpha=1)+
  scale_color_manual("", values=cpList$smoothed.bright)
p2 <- plotReducedDim(spe_pseudo, dimred="PCA_1663", ncomponents=2, colour_by = "condition", shape_by= "sex", point_alpha=1)+
  scale_color_manual("", values=cpList$dx.pal)+
  scale_shape_manual(values=c(21,23))
b1 = ggplot_build(p1)
b1$data[[1]]$size = .5
b1$data[[1]]$shape = 21
# because of overlap of points I like empty fill better
b1$data[[1]]$fill = NA

b2 = ggplot_build(p2)
b2$data[[1]]$size = .5
b2$data[[1]]$fill = NA
## alt to fill with lighter version of dx.pal color
## based on the small size of the points I actually think i like empty better
#b2$data[[1]]$fill = as.character(factor(b2$data[[1]]$colour, levels=c(cpList$dx.pal[[1]], cpList$dx.pal[[2]], cpList$dx.pal[[3]]),
#                             labels=c(colorRampPalette(c("white",cpList$dx.pal[[1]]))(20)[[12]],
#                                      colorRampPalette(c("white",cpList$dx.pal[[2]]))(20)[[12]],
#                                      colorRampPalette(c("white",cpList$dx.pal[[3]]))(20)[[12]])))

ggsave(file="plots/publication/Figure1/PC1-PC2_UMAP.pdf",
       grid.arrange(ggplot_gtable(b1), ggplot_gtable(b2), ncol=1),
       width=3, height=4)

# set up for different variable sets
colnames(colData(spe_pseudo))[grep("smoothed", colnames(colData(spe_pseudo)))] = "cluster"
colnames(colData(spe_pseudo))[grep("subsets_mito_percent", colnames(colData(spe_pseudo)))] = "chrM_ratio"
spe_pseudo$pc3 = reducedDim(spe_pseudo, "PCA_1663")[,3]

cdata = read.csv("raw-data/sample_info/Page_BMI_Smoking_CrossDisorder_120825.csv")
cdata = filter(cdata, BrNum %in% spe_pseudo$brnum)
sdata = distinct(as.data.frame(colData(spe_pseudo)[,c("brnum","sample_id","condition","sex","age","PMI","RIN")]))
sdata = left_join(sdata, cdata[,c("BrNum","BMI","Smoking")], by=c("brnum"="BrNum"))

# update colData
new.cdata = merge(colData(spe_pseudo), sdata, sort=F)
stopifnot(identical(spe_pseudo$total, new.cdata$total))
colData(spe_pseudo) <- new.cdata

# pca with new variables
bio.vars = c("cluster","sum","detected","nspots","chrM_ratio")
bio.colors = RColorBrewer::brewer.pal(n=length(bio.vars), "Set1")
names(bio.colors) = bio.vars

exp.vars = c("sample_id", "condition", "sex", "slide", "seq","pc3")
exp.colors = RColorBrewer::brewer.pal(n=length(exp.vars), "Dark2")
names(exp.colors) = exp.vars
exp.vars = exp.vars[1:(length(exp.vars)-1)]
exp.colors = exp.colors[1:(length(exp.colors)-1)]

donor.vars = c("age", "BMI", "Smoking", "RIN", "PMI")
donor.colors = c("#E4775D","#A0C255","#EEBC4A","grey","#93D3F6")
names(donor.colors) = donor.vars


p1 <- plotExplanatoryPCs(spe_pseudo, dimred="PCA_1663", npcs_to_plot=10, variables=exp.vars)+
        scale_y_continuous()+
	scale_color_manual("", values= exp.colors)+
	scale_x_continuous(breaks=c(2,4,6,8,10))+
	labs(subtitle="Experimental variables", y="% PC variance explained")+
	theme(legend.position="bottom")

p2 <- plotExplanatoryPCs(spe_pseudo, dimred="PCA_1663", npcs_to_plot=10, variables=bio.vars)+
        scale_y_continuous()+
        scale_color_manual("", values=bio.colors)+
	scale_x_continuous(breaks=c(2,4,6,8,10))+
        labs(subtitle="Domain variables", y="% PC variance explained")+
	theme(legend.position="bottom")

p3 <- plotExplanatoryPCs(spe_pseudo, dimred="PCA_1663", npcs_to_plot=10, variables=donor.vars)+
        scale_y_continuous()+
        scale_color_manual("", values=donor.colors)+
	scale_x_continuous(breaks=c(2,4,6,8,10))+
        labs(subtitle="Donor variables", y="% PC variance explained")+
	theme(legend.position="bottom")

ggsave(file="plots/publication/Figure1/supp_pca-variance-explained.pdf",
	grid.arrange(p2, p1, p3, ncol=3), height=2.5, width=8)


cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()

