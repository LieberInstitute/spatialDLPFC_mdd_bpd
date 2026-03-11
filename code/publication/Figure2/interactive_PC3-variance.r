library(SpatialExperiment)
library(scater)
library(dplyr)
library(ggplot2)
set.seed(123)

load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
colnames(colData(spe_pseudo))[grep("subsets_mito_percent", colnames(colData(spe_pseudo)))] = "chrM_ratio"
spe_pseudo$pc3 = reducedDim(spe_pseudo, "PCA_1663")[,3]
spe_pseudo$domain = spe_pseudo$smoothed_k9_1663

#add extra covars
cdata = read.csv("raw-data/sample_info/Page_BMI_Smoking_CrossDisorder_120825.csv")
sdata = distinct(as.data.frame(colData(spe_pseudo)[,c("brnum","sample_id","condition","sex","age","PMI","RIN")]))
sdata = left_join(sdata, cdata[,c("BrNum","BMI","Smoking")], by=c("brnum"="BrNum"))
new.cdata = merge(colData(spe_pseudo), sdata, sort=F)
stopifnot(identical(spe_pseudo$total, new.cdata$total))
colData(spe_pseudo) <- new.cdata

colnames(colData(spe_pseudo))

# pca with new variables
bio.vars = c("domain","sum","detected","nspots","chrM_ratio")
bio.colors = RColorBrewer::brewer.pal(n=length(bio.vars), "Set1")
names(bio.colors) = bio.vars

exp.vars = c("sample_id", "condition", "sex", "slide", "seq","pc3")
exp.colors = RColorBrewer::brewer.pal(n=length(exp.vars), "Dark2")
names(exp.colors) = exp.vars

donor.vars = c("age", "BMI", "Smoking", "RIN")#, "PMI")
donor.colors = c("#E4775D","#A0C255","#EEBC4A","grey")#,"#93D3F6")
names(donor.colors) = donor.vars



var.m <- getExplanatoryPCs(spe_pseudo, dimred="PCA_1663", variables= c(bio.vars, exp.vars, donor.vars))

tmp2 = rbind(cbind.data.frame("PC"=rep("PC1", ncol(var.m)), variance=var.m["PC1",]),
             cbind.data.frame("PC"=rep("PC2", ncol(var.m)), variance=var.m["PC2",]),
             cbind.data.frame("PC"=rep("PC3", ncol(var.m)), variance=var.m["PC3",]))
tmp2$variable = rep(colnames(var.m), 3)

tmp2$variable_f = factor(tmp2$variable, levels=rev(c(bio.vars, exp.vars, donor.vars)))

p1 <- ggplot(filter(tmp2, variable_f!="pc3"), aes(y=variable_f, x=variance, fill=variable))+
  geom_bar(stat="identity", position="stack")+
  scale_fill_manual(values=c(bio.colors, exp.colors, donor.colors), guide="none")+
  facet_wrap(vars(PC), ncol=3)+
  labs(x="variance explained", y="experimental variables", title="PC3: continuous batch effect")+
  theme_minimal()+theme(plot.margin = margin(.5,1,.5,.5, "cm"), 
                        axis.text.x=element_text(size=8),
                        panel.spacing = unit(.5, "cm"))

ggsave(file="plots/publication/Figure2/PC3_variance.png", 
       p1, bg="white",
       width=5, height=5)

