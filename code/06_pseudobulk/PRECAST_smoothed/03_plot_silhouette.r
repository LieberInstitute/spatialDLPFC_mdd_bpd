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

#PRECAST smoothed
load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo-with-lowUMI_sample-smoothed-n1663-k9_norm.Rdata")

table(spe_pseudo$smoothed_k9_1663, useNA="ifany")

spe_pseudo$cond_sex = factor(paste(spe_pseudo$condition, spe_pseudo$sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"))
table(as.data.frame(colData(spe_pseudo))[,c("cond_sex","smoothed_k9_1663")], useNA="ifany")

color.palette = c(cpList$smoothed.bright, "low UMI"="grey50")
sil.results <- as.data.frame(approxSilhouette(reducedDim(spe_pseudo, "PCA_1663"), clusters=spe_pseudo$smoothed_k9_1663))
sil.results$closest <- factor(ifelse(sil.results$width > 0, as.character(sil.results$cluster), as.character(sil.results$other)))
sil.results$closest <- factor(sil.results$closest, levels=levels(spe_pseudo$smoothed_k9_1663))

p1 <- ggplot(sil.results, aes(x=cluster, y=width, colour=closest))+
  ggbeeswarm::geom_quasirandom()+scale_color_manual(values=color.palette)+
  labs(x="PRECAST (smoothed) cluster", title="PRECAST (smoothed)", subtitle="with low UMI cluster")+
  theme_minimal()+theme(text=element_text(size=10))


#PRECAST smoothed
load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm.Rdata")

table(spe_pseudo$smoothed_k9_1663, useNA="ifany")

spe_pseudo$cond_sex = factor(paste(spe_pseudo$condition, spe_pseudo$sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"))
table(as.data.frame(colData(spe_pseudo))[,c("cond_sex","smoothed_k9_1663")], useNA="ifany")


sil.results <- as.data.frame(approxSilhouette(reducedDim(spe_pseudo, "PCA_1663"), clusters=spe_pseudo$smoothed_k9_1663))
sil.results$closest <- factor(ifelse(sil.results$width > 0, as.character(sil.results$cluster), as.character(sil.results$other)))
sil.results$closest <- factor(sil.results$closest, levels=levels(spe_pseudo$smoothed_k9_1663))

p2 <- ggplot(sil.results, aes(x=cluster, y=width, colour=closest))+
  ggbeeswarm::geom_quasirandom()+scale_color_manual(values=cpList$smoothed.bright)+
  labs(x="PRECAST (smoothed) cluster", title="PRECAST (smoothed)", subtitle="without low UMI cluster, unfiltered")+
  theme_minimal()+theme(text=element_text(size=10))


#PRECAST smoothed filtered
load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")

table(spe_pseudo$smoothed_k9_1663, useNA="ifany")

spe_pseudo$cond_sex = factor(paste(spe_pseudo$condition, spe_pseudo$sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"))
table(as.data.frame(colData(spe_pseudo))[,c("cond_sex","smoothed_k9_1663")], useNA="ifany")


sil.results <- as.data.frame(approxSilhouette(reducedDim(spe_pseudo, "PCA_1663"), clusters=spe_pseudo$smoothed_k9_1663))
sil.results$closest <- factor(ifelse(sil.results$width > 0, as.character(sil.results$cluster), as.character(sil.results$other)))
sil.results$closest <- factor(sil.results$closest, levels=levels(spe_pseudo$smoothed_k9_1663))

p3 <- ggplot(sil.results, aes(x=cluster, y=width, colour=closest))+
  ggbeeswarm::geom_quasirandom()+scale_color_manual(values=cpList$smoothed.bright)+
  labs(x="PRECAST (smoothed) cluster", title="PRECAST (smoothed)", subtitle="without low UMI cluster, filtered")+
  theme_minimal()+theme(text=element_text(size=10))


ggsave(file="plots/06_pseudobulk/PRECAST_smoothed/PRECAST_silhouettes.png", 
       grid.arrange(p1, p2, p3, layout_matrix=rbind(c(1,NA),c(2,3))),
       bg="white", height=8, width=8)
cat("\nSaved silhouette plots to: plots/06_pseudobulk/PRECAST_smoothed/PRECAST_silhouettes.png\n")

cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
