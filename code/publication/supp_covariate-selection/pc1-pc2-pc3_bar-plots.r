setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(pheatmap)
	library(scater)
	library(dplyr)
	library(ggplot2)
})


load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
colnames(colData(spe_pseudo))[grep("subsets_mito_percent", colnames(colData(spe_pseudo)))] = "chrM_ratio"
spe_pseudo$pc3 = reducedDim(spe_pseudo, "PCA_1663")[,3]
spe_sm <- spe_pseudo

#need to remove age from colData because of significant digits change
demo = read.csv("processed-data/publication/demographics.csv")
new.cdata = merge(colData(spe_sm)[,setdiff(colnames(colData(spe_sm)), c("age","RIN"))], demo, sort=F)
stopifnot(identical(spe_sm$total, new.cdata$total))
colData(spe_sm) <- new.cdata




load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
colnames(colData(spe_pseudo))[grep("subsets_mito_percent", colnames(colData(spe_pseudo)))] = "chrM_ratio"
spe_pseudo$pc3 = reducedDim(spe_pseudo, "PCA_1663")[,3]
spe_se <- spe_pseudo
new.cdata = merge(colData(spe_se)[,setdiff(colnames(colData(spe_se)), c("age","RIN"))], demo, sort=F)
stopifnot(identical(spe_se$total, new.cdata$total))
colData(spe_se) <- new.cdata


# define variable groups
core.vars = c("condition", "sex", "sample_id")
other.vars = c("nspots", "chrM_ratio", "pc3","age", "BMI", "RIN", "Smoking", "slide", "seq")


# PC variance explained
var.pcs <- getExplanatoryPCs(spe_sm, dimred="PCA_1663", variables= c("smoothed_k9_1663", core.vars, 
                                                                     "detected","sum", other.vars))
tmp = rbind(cbind.data.frame("PC"=rep("PC1", ncol(var.pcs)), variance=var.pcs["PC1",]),
            cbind.data.frame("PC"=rep("PC2", ncol(var.pcs)), variance=var.pcs["PC2",]),
            cbind.data.frame("PC"=rep("PC3", ncol(var.pcs)), variance=var.pcs["PC3",]))
tmp$variable = rep(colnames(var.pcs), 3)
tmp$variable_f = factor(tmp$variable, levels=rev(c("smoothed_k9_1663","detected","nspots","sum","chrM_ratio","sample_id","slide","seq",
                                                   "age","Smoking","RIN","condition","sex","BMI","pc3")),
                        labels=rev(c("domain","detected","nspots","sum","chrM_ratio","sample_id","slide","seq",
                                     "age","Smoking","RIN","condition","sex","BMI","pc3")))


p1 <- ggplot(filter(tmp, variable_f!="pc3"), aes(y=variable_f, x=variance))+
  geom_bar(stat="identity")+
  #scale_fill_manual(values=c(bio.colors, exp.colors, donor.colors), guide="none")+
  facet_wrap(vars(PC), ncol=3)+
  coord_cartesian(xlim=c(0,100))+
  labs(x="variance explained", y="experimental variables", title="PC3: PRECAST (smoothed)")+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(),
                        axis.text.x=element_text(size=8), plot.margin = margin(.5,1,.5,.5, "cm"), 
                        panel.spacing = unit(.5, "cm"))


var.pcs2 <- getExplanatoryPCs(spe_se, dimred="PCA_1663", variables= c("seurat_label", core.vars, 
                                                                      "detected","sum", other.vars))
tmp2 = rbind(cbind.data.frame("PC"=rep("PC1", ncol(var.pcs2)), variance=var.pcs2["PC1",]),
             cbind.data.frame("PC"=rep("PC2", ncol(var.pcs2)), variance=var.pcs2["PC2",]),
             cbind.data.frame("PC"=rep("PC3", ncol(var.pcs2)), variance=var.pcs2["PC3",]))
tmp2$variable = rep(colnames(var.pcs2), 3)
tmp2$variable_f = factor(tmp2$variable, levels=rev(c("seurat_label","detected","nspots","sum","chrM_ratio","sample_id","slide","seq",
                                                     "age","Smoking","RIN","condition","sex","BMI","pc3")),
                         labels=rev(c("domain","detected","nspots","sum","chrM_ratio","sample_id","slide","seq",
                                      "age","Smoking","RIN","condition","sex","BMI","pc3")))


p2 <- ggplot(filter(tmp2, variable_f!="pc3"), aes(y=variable_f, x=variance))+
  geom_bar(stat="identity")+
  facet_wrap(vars(PC), ncol=3)+
  coord_cartesian(xlim=c(0,100))+
  labs(x="variance explained", y="experimental variables", title="PC3: Seurat label")+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(),
                        axis.text.x=element_text(size=8), plot.margin = margin(.5,1,.5,.5, "cm"), 
                        panel.spacing = unit(.5, "cm"))

pdf(file="plots/publication/supp_covariate-selection/top3-pc3_bar-plots.pdf")
p1
p2
dev.off()

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()

