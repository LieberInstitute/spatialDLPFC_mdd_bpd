setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(bluster)
})

set.seed(123)

cpList <- readRDS("plots/colorPalettes.rds")
names(cpList$transfer.bright)[1] = "M.V"
names(cpList$low.res.bright)[1] = "M.V"
fill.palette = c('M.V'="#911223", 'Ast.L1'="#cfa45c", 'Ast.Nrn'= "#F5D29E",
                 'L2'= "#5D9940", 'L3'= "#5095CD",'L4'= "#c2cfcf",'Inhb'= "#9377AC",
                 'L5'= "#ddc94e", 'L6'= "#E45C5F", 'WM'= "#D1C4B0")

cat("\nCustom cluster pseudobulk, before filtering...\n")
load("processed-data/06_pseudobulk/custom_cluster/spe_n119_pseudo_sample-custom-cluster_norm.Rdata")
spe_pseudo$custom_cluster = factor(spe_pseudo$custom_cluster, levels=c("Micro.Vasc","Astro.L1","Astro.Nrn","L2","L3","L4","Inhb","L5","L6","WM"),
                                   labels=c("M.V","Ast.L1","Ast.Nrn","L2","L3","L4","Inhb","L5","L6","WM"))

table(spe_pseudo$custom_cluster, useNA="ifany")

spe_pseudo$cond_sex = factor(paste(spe_pseudo$condition, spe_pseudo$sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"))
table(as.data.frame(colData(spe_pseudo))[,c("cond_sex","custom_cluster")], useNA="ifany")

sil.results <- as.data.frame(approxSilhouette(reducedDim(spe_pseudo, "PCA_1663"), clusters=spe_pseudo$custom_cluster))
sil.results$closest <- factor(ifelse(sil.results$width > 0, as.character(sil.results$cluster), as.character(sil.results$other)))
sil.results$closest <- factor(sil.results$closest, levels=levels(spe_pseudo$custom_cluster))

p1 <- ggplot(sil.results, aes(x=cluster, y=width, colour=closest))+
  ggbeeswarm::geom_quasirandom()+scale_color_manual(values=fill.palette)+
  labs(x="cluster", title="Custom clusters", subtitle="before QC filter")+
  theme_minimal()+theme(text=element_text(size=10))

#QC filtered
cat("\nCustom cluster pseudobulk, after filtering...\n")
load("processed-data/06_pseudobulk/custom_cluster/spe_n119_pseudo_sample-custom-cluster_norm-filt.Rdata")
spe_pseudo$custom_cluster = factor(spe_pseudo$custom_cluster, levels=c("Micro.Vasc","Astro.L1","Astro.Nrn","L2","L3","L4","Inhb","L5","L6","WM"),
                                   labels=c("M.V","Ast.L1","Ast.Nrn","L2","L3","L4","Inhb","L5","L6","WM"))

table(spe_pseudo$custom_cluster, useNA="ifany")

spe_pseudo$cond_sex = factor(paste(spe_pseudo$condition, spe_pseudo$sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"))
table(as.data.frame(colData(spe_pseudo))[,c("cond_sex","custom_cluster")], useNA="ifany")

sil.results <- as.data.frame(approxSilhouette(reducedDim(spe_pseudo, "PCA_1663"), clusters=spe_pseudo$custom_cluster))
sil.results$closest <- factor(ifelse(sil.results$width > 0, as.character(sil.results$cluster), as.character(sil.results$other)))
sil.results$closest <- factor(sil.results$closest, levels=levels(spe_pseudo$custom_cluster))

p2 <- ggplot(sil.results, aes(x=cluster, y=width, colour=closest))+
  ggbeeswarm::geom_quasirandom()+scale_color_manual(values=fill.palette)+
  labs(x="cluster", title="Custom clusters", subtitle="after QC filter")+
  theme_minimal()+theme(text=element_text(size=10))


ggsave(file="plots/06_pseudobulk/custom_cluster/custom-cluster_silhouettes.png", 
       grid.arrange(p1, p2, ncol=2),
       bg="white", height=4, width=8)
cat("\nSaved silhouette plots to: plots/06_pseudobulk/custom_cluster/custom-cluster_silhouettes.png\n")


cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
