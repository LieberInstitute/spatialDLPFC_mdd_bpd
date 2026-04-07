setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(scater)
  library(gridExtra)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(bluster)
})

set.seed(123)


cpList = readRDS("plots/colorPalettes.rds")

load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
dim(spe_pseudo)
min(spe_pseudo$detected)

p1 <- plotReducedDim(spe_pseudo, dimred="PCA_1663", ncomponents=2, colour_by = "seurat_label", point_alpha=1)+
  scale_color_manual("", values=cpList$transfer.bright, labels=c("M.V","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))
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

ggsave(file="plots/publication/Figure3/supp_PC1-PC2_UMAP.pdf",
       grid.arrange(ggplot_gtable(b1), ggplot_gtable(b2), ncol=1),
       width=3, height=4)


# silhouette
sil.results3 <- as.data.frame(approxSilhouette(reducedDim(spe_pseudo, "PCA_1663"), clusters=spe_pseudo$seurat_label))
sil.results3$closest <- factor(ifelse(sil.results3$width > 0, as.character(sil.results3$cluster), as.character(sil.results3$other)))
sil.results3$closest <- factor(sil.results3$closest, levels=levels(spe_pseudo$seurat_label))

p3 <- ggplot(sil.results3, aes(x=cluster, y=width, colour=closest))+
  ggbeeswarm::geom_quasirandom()+geom_hline(aes(yintercept=0), lty=3, linewidth=1)+
  scale_color_manual("closest\nseurat_label", values=cpList$transfer.bright)+
  scale_x_discrete(labels=c("Micro/\nVasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))+
  labs(x="", y="silhouette width")+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank())

# pseudobulk qc of unfiltered
spe_save = spe_pseudo
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm.Rdata")

p4 <- ggplot(as.data.frame(colData(spe_pseudo)), 
             aes(x=seurat_label, y=detected, color=seurat_label))+
  ggbeeswarm::geom_quasirandom()+
  scale_color_manual(values=cpList$transfer.bright)+
  labs(x="", y="pseudobulk detected genes")+
  scale_x_discrete(labels=c("Micro/\nVasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))+
  geom_hline(aes(yintercept=10000), lty=2, linewidth=1)+
  geom_hline(aes(yintercept=8000), lty=2, linewidth=1)+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank())


ggsave(file="plots/publication/Figure3/supp_seurat-pseudobulk.pdf", 
	grid.arrange(p4, p3, ncol=1), height=5, width=5)



# set up for different variable sets
spe_pseudo = spe_save
colnames(colData(spe_pseudo))[grep("seurat_label", colnames(colData(spe_pseudo)))] = "cluster"
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

ggsave(file="plots/publication/Figure3/supp_seurat-label_pca-variance-explained.pdf",
        grid.arrange(p2, p1, p3, ncol=3), height=2.5, width=8)


## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
