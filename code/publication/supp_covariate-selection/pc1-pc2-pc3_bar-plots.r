setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(pheatmap)
	library(scater)
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
})
set.seed(123)

load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
colnames(colData(spe_pseudo))[grep("subsets_mito_percent", colnames(colData(spe_pseudo)))] = "chrM_ratio"
spe_pseudo$pc3 = reducedDim(spe_pseudo, "PCA_1663")[,3]
spe_sm <- spe_pseudo

#add additional covars
demo = read.csv("processed-data/publication/supp_tables/demographics.csv")

new.cdata = merge(colData(spe_sm)[,c("sample_id","brnum","condition","sex","smoothed_k9_1663","nspots","chrM_ratio","pc3","detected","sum","slide","seq")],
        demo[,c("sample_id","brnum","sex","age","RIN","BMI","Smoking")], sort=F)
stopifnot(identical(spe_sm$pc3, new.cdata$pc3))
colData(spe_sm) <- new.cdata




load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
colnames(colData(spe_pseudo))[grep("subsets_mito_percent", colnames(colData(spe_pseudo)))] = "chrM_ratio"
spe_pseudo$pc3 = reducedDim(spe_pseudo, "PCA_1663")[,3]
spe_se <- spe_pseudo

new.cdata = merge(colData(spe_se)[,c("sample_id","brnum","condition","sex","seurat_label","nspots","chrM_ratio","pc3","detected","sum","slide","seq")],
        demo[,c("sample_id","brnum","sex","age","RIN","BMI","Smoking")], sort=F)
stopifnot(identical(spe_se$pc3, new.cdata$pc3))
colData(spe_se) <- new.cdata


# define variable groups
core.vars = c("condition", "sex", "sample_id")
other.vars = c("nspots", "chrM_ratio","age", "BMI", "RIN", "Smoking", "slide", "seq")


# PC variance explained
var.pcs <- getExplanatoryPCs(spe_sm, dimred="PCA_1663", variables= c("smoothed_k9_1663", core.vars, 
                                                                     "detected","sum", other.vars))
tmp = rbind(cbind.data.frame("PC"=rep("PC1", ncol(var.pcs)), variance=var.pcs["PC1",]),
            cbind.data.frame("PC"=rep("PC2", ncol(var.pcs)), variance=var.pcs["PC2",]),
            cbind.data.frame("PC"=rep("PC3", ncol(var.pcs)), variance=var.pcs["PC3",]))
tmp$variable = rep(colnames(var.pcs), 3)
tmp$variable_f = factor(tmp$variable, levels=rev(c("smoothed_k9_1663","detected","nspots","sum","chrM_ratio","sample_id","slide","seq",
                                                   "age","Smoking","RIN","condition","sex","BMI")),
                        labels=rev(c("domain","detected","nspots","sum","chrM_ratio","sample_id","slide","seq",
                                     "age","Smoking","RIN","condition","sex","BMI")))

tmp$PC = factor(tmp$PC, levels=c("PC3","PC2","PC1"))

p1 <- ggplot(tmp, aes(y=variable_f, x=variance, fill=PC))+
  geom_bar(stat="identity", position="dodge", width=.9)+
  scale_fill_manual(values=c("#e41a1c","#377eb8","#a65628"))+
  coord_cartesian(xlim=c(0,100))+
  labs(x="variance explained", y="experimental variables", title="domain-SP")+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(),
                        text=element_text(size=6), legend.key.size = unit(6,"pt"))

# now domain-CT
var.pcs2 <- getExplanatoryPCs(spe_se, dimred="PCA_1663", variables= c("seurat_label", core.vars, 
                                                                      "detected","sum", other.vars))
tmp2 = rbind(cbind.data.frame("PC"=rep("PC1", ncol(var.pcs2)), variance=var.pcs2["PC1",]),
             cbind.data.frame("PC"=rep("PC2", ncol(var.pcs2)), variance=var.pcs2["PC2",]),
             cbind.data.frame("PC"=rep("PC3", ncol(var.pcs2)), variance=var.pcs2["PC3",]))
tmp2$variable = rep(colnames(var.pcs2), 3)
tmp2$variable_f = factor(tmp2$variable, levels=rev(c("seurat_label","detected","nspots","sum","chrM_ratio","sample_id","slide","seq",
                                                     "age","Smoking","RIN","condition","sex","BMI")),
                         labels=rev(c("domain","detected","nspots","sum","chrM_ratio","sample_id","slide","seq",
                                      "age","Smoking","RIN","condition","sex","BMI")))

tmp2$PC = factor(tmp2$PC, levels=c("PC3","PC2","PC1"))

p2 <- ggplot(tmp2, aes(y=variable_f, x=variance, fill=PC))+
  geom_bar(stat="identity", position="dodge", width=.9)+
  scale_fill_manual(values=c("#e41a1c","#377eb8","#a65628"))+
  coord_cartesian(xlim=c(0,100))+
  labs(x="variance explained", y="experimental variables", title="domain-CT")+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(),
                        text=element_text(size=6), legend.key.size = unit(6,"pt"))

pdf(file="plots/publication/supp_covariate-selection/top3-pcs_bar-plots.pdf", height=3, width=3)
grid.arrange(p1+theme(legend.position="none"), p2+theme(legend.position="none"), ncol=2)
grid.arrange(p1, p2, ncol=2)
dev.off()

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()

