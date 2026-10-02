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


ggsave(file="plots/publication/supp_clustering_Seurat/seurat-label_pseudobulk-qc_violin.pdf", 
	grid.arrange(p4, p3, ncol=1), height=5, width=5)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
