library(SpatialExperiment)
library(scater)
library(dplyr)
library(ggplot2)
set.seed(123)

load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-no-lowUMI_sample-seurat-pc30_norm-filt.Rdata")

exp.vars = c("seurat_label","sample_id",
             "slide","seq",
             "condition","sex",
             "age","PMI","RIN",
             "sum","detected","nspots")
exp.vars.colors = c("#FB8072", "#80B1D3",
                    "#BC80BD", "#B3DE69",
                    "#8DD3C7", "#FDB462",
                    "#FCCDE5", "#BEBADA", "#FFFFB3",
                    "#CCEBC5", "#D9D9D9","black")
names(exp.vars.colors) = exp.vars
var.m <- getExplanatoryPCs(spe_pseudo, dimred="PCA_1663", variables= exp.vars)

head(var.m)

tmp = cbind.data.frame("PC"=rep("PC3", ncol(var.m)), var.m["PC3",])
colnames(tmp) = c("PC","variance")
tmp$variable = rownames(tmp)


tmp[tmp$variable=="seurat_label","variable"] = "cluster"
names(exp.vars.colors)[1] = "cluster"
tmp$variable_f = factor(tmp$variable, levels=tmp$variable[order(tmp$variance)])

p1 <- ggplot(tmp, aes(y=variable_f, x=variance, fill=variable))+
  geom_bar(stat="identity", position="stack")+
  scale_fill_manual(values=exp.vars.colors, guide="none")+
  labs(x="variance explained", y="experimental variables", title="PC3")+
  theme_minimal()+theme(plot.margin = margin(.5,1,.5,.5, "cm"))

ggsave(file="plots/publication/PC3_variance.png", 
       p1, bg="white",
       width=5, height=5)

tmp2 = rbind(cbind.data.frame("PC"=rep("PC1", ncol(var.m)), variance=var.m["PC1",]),
             cbind.data.frame("PC"=rep("PC2", ncol(var.m)), variance=var.m["PC2",]),
             cbind.data.frame("PC"=rep("PC3", ncol(var.m)), variance=var.m["PC3",]))
tmp2$variable = rep(colnames(var.m), 3)

tmp2[tmp2$variable=="seurat_label","variable"] = "cluster"
names(exp.vars.colors)[1] = "cluster"
tmp2$variable_f = factor(tmp2$variable, levels=tmp$variable[order(tmp$variance)])

p1 <- ggplot(tmp2, aes(y=variable_f, x=variance, fill=variable))+
  geom_bar(stat="identity", position="stack")+
  scale_fill_manual(values=exp.vars.colors, guide="none")+
  facet_wrap(vars(PC), ncol=3)+
  labs(x="variance explained", y="experimental variables", title="PC3: continuous batch effect")+
  theme_minimal()+theme(plot.margin = margin(.5,1,.5,.5, "cm"), 
                        axis.text.x=element_text(size=8),
                        panel.spacing = unit(.5, "cm"))

ggsave(file="plots/publication/PC3_variance.png", 
       p1, bg="white",
       width=5, height=5)
