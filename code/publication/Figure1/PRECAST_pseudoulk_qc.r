setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(dplyr)
  library(ggplot2)
  library(bluster)
  library(gridExtra)
})


set.seed(123)

cpList = readRDS("plots/colorPalettes.rds")

#load spe
load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo-with-lowUMI_sample-smoothed-n1663-k9_norm.Rdata")

# umi counts with low UMI cluster
p1 <- ggplot(as.data.frame(colData(spe_pseudo)), 
             aes(x=smoothed_k9_1663, y=sum, color=smoothed_k9_1663))+
  ggbeeswarm::geom_quasirandom()+scale_y_log10()+
  scale_color_manual(values=c("low UMI"="grey50", cpList$smoothed.bright))+
  labs(title="QC", x="", y="pseudobulk total UMI")+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank())

# silhouette plot with low UMI cluster
sil.results <- as.data.frame(approxSilhouette(reducedDim(spe_pseudo, "PCA_1663"), clusters=spe_pseudo$smoothed_k9_1663))
sil.results$closest <- factor(ifelse(sil.results$width > 0, as.character(sil.results$cluster), as.character(sil.results$other)))
sil.results$closest <- factor(sil.results$closest, levels=levels(spe_pseudo$smoothed_k9_1663))

p2 <- ggplot(sil.results, aes(x=cluster, y=width, colour=closest))+
  ggbeeswarm::geom_quasirandom()+geom_hline(aes(yintercept=0), lty=3, linewidth=1)+
  scale_color_manual("closest\nsmoothed_k9_1663", values=c("low UMI"="grey50", cpList$smoothed.bright))+
  labs(title="Silhouette", x="", y="silhouette width")+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank())


# no low UMI cluster
load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm.Rdata")

# detected genes counts without low UMI cluster
p3 <- ggplot(as.data.frame(colData(spe_pseudo)), 
             aes(x=smoothed_k9_1663, y=detected, color=smoothed_k9_1663))+
  ggbeeswarm::geom_quasirandom()+
  scale_color_manual(values=cpList$smoothed.bright)+
  labs(title="QC", x="", y="pseudobulk detected genes")+
  geom_hline(aes(yintercept=7500), lty=2, linewidth=1)+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank())

# silhouette plot with low UMI cluster
sil.results2 <- as.data.frame(approxSilhouette(reducedDim(spe_pseudo, "PCA_1663"), clusters=spe_pseudo$smoothed_k9_1663))
sil.results2$closest <- factor(ifelse(sil.results2$width > 0, as.character(sil.results2$cluster), as.character(sil.results2$other)))
sil.results2$closest <- factor(sil.results2$closest, levels=levels(spe_pseudo$smoothed_k9_1663))

p4 <- ggplot(sil.results2, aes(x=cluster, y=width, colour=closest))+
  ggbeeswarm::geom_quasirandom()+geom_hline(aes(yintercept=0), lty=3, linewidth=1)+
  scale_color_manual("closest\nsmoothed_k9_1663", values=cpList$smoothed.bright)+
  labs(title="Pre-filter", x="", y="silhouette width")+
  coord_cartesian(ylim=c(-.75, .5))+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
                        legend.position="none")

# no low UMI cluster, filtered
load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
sil.results3 <- as.data.frame(approxSilhouette(reducedDim(spe_pseudo, "PCA_1663"), clusters=spe_pseudo$smoothed_k9_1663))
sil.results3$closest <- factor(ifelse(sil.results3$width > 0, as.character(sil.results3$cluster), as.character(sil.results3$other)))
sil.results3$closest <- factor(sil.results3$closest, levels=levels(spe_pseudo$smoothed_k9_1663))

p5 <- ggplot(sil.results3, aes(x=cluster, y=width, colour=closest))+
  ggbeeswarm::geom_quasirandom()+geom_hline(aes(yintercept=0), lty=3, linewidth=1)+
  scale_color_manual("closest\nsmoothed_k9_1663", values=cpList$smoothed.bright)+
  labs(title="Post-filter", x="", y="silhouette width")+
  coord_cartesian(ylim=c(-.75, .5))+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
                        legend.position="none")


pdf(file="plots/publication/Figure1/supp_PRECAST-pseudobulk.pdf", height=5, width=5)
grid.arrange(p1, p2, ncol=1)
grid.arrange(p3, p4, p5, layout_matrix=rbind(c(1,1),c(2,3)))
dev.off()


## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
