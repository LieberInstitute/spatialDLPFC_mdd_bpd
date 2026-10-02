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
### no detected genes filter at this point
table(spe_pseudo$smoothed_k9_1663, useNA="ifany")
#.L1      L2    L3.4      L5      L6      WM low UMI 
#119     119     119     119     119     105     110 
spe_pseudo$cond_sex = factor(paste(spe_pseudo$condition, spe_pseudo$sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"))
table(as.data.frame(colData(spe_pseudo))[,c("cond_sex","smoothed_k9_1663")], useNA="ifany")
#.........smoothed_k9_1663
#cond_sex L1 L2 L3.4 L5 L6 WM low UMI
#...NTC F 20 20   20 20 20 16      18
#...NTC M 20 20   20 20 20 17      20
#...MDD F 20 20   20 20 20 16      18
#...MDD M 19 19   19 19 19 18      18
#...BPD F 20 20   20 20 20 19      17
#...BPD M 20 20   20 20 20 19      19


color.palette = c(cpList$smoothed.bright, "low UMI"="grey50")
sil.results <- as.data.frame(approxSilhouette(reducedDim(spe_pseudo, "PCA_1663"), clusters=spe_pseudo$smoothed_k9_1663))
sil.results$closest <- factor(ifelse(sil.results$width > 0, as.character(sil.results$cluster), as.character(sil.results$other)))
sil.results$closest <- factor(sil.results$closest, levels=levels(spe_pseudo$smoothed_k9_1663))

p1 <- ggplot(sil.results, aes(x=cluster, y=width, colour=closest))+
  ggbeeswarm::geom_quasirandom()+scale_color_manual(values=color.palette)+
  labs(x="PRECAST (smoothed) cluster", title="PRECAST (smoothed)", subtitle="with low UMI cluster")+
  theme_minimal()+theme(text=element_text(size=10))

spe_sm_low.umi = spe_pseudo

#PRECAST smoothed
load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
### detected genes filter is 10k across the board
table(spe_pseudo$smoothed_k9_1663, useNA="ifany")
#.L1   L2 L3.4   L5   L6   WM 
#117  119  119  119  118   93 
spe_pseudo$cond_sex = factor(paste(spe_pseudo$condition, spe_pseudo$sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"))
table(as.data.frame(colData(spe_pseudo))[,c("cond_sex","smoothed_k9_1663")], useNA="ifany")
#.........smoothed_k9_1663
#cond_sex L1 L2 L3.4 L5 L6 WM
#...NTC F 20 20   20 20 20 14
#...NTC M 18 20   20 20 20 16
#...MDD F 20 20   20 20 19 15
#...MDD M 19 19   19 19 19 15
#...BPD F 20 20   20 20 20 14
#...BPD M 20 20   20 20 20 19

sil.results <- as.data.frame(approxSilhouette(reducedDim(spe_pseudo, "PCA_1663"), clusters=spe_pseudo$smoothed_k9_1663))
sil.results$closest <- factor(ifelse(sil.results$width > 0, as.character(sil.results$cluster), as.character(sil.results$other)))
sil.results$closest <- factor(sil.results$closest, levels=levels(spe_pseudo$smoothed_k9_1663))

p2 <- ggplot(sil.results, aes(x=cluster, y=width, colour=closest))+
  ggbeeswarm::geom_quasirandom()+scale_color_manual(values=cpList$smoothed.bright)+
  labs(x="PRECAST (smoothed) cluster", title="PRECAST (smoothed)", subtitle="without low UMI cluster")+
  theme_minimal()+theme(text=element_text(size=10))

spe_sm = spe_pseudo


#label transfer
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
### detected genes filter is 10k across the board
spe_pseudo$seurat_label_f = factor(as.character(spe_pseudo$seurat_label), 
                                   levels=c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"),
                                   labels=c("M.V","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))
table(spe_pseudo$seurat_label_f, useNA="ifany")
#Micro.Vasc      Astro       L2.3         L4       Inhb         L5         L6 
#114        119        119        119        119        119        119 
#Oligo 
#105 
spe_pseudo$cond_sex = factor(paste(spe_pseudo$condition, spe_pseudo$sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"))
table(as.data.frame(colData(spe_pseudo))[,c("cond_sex","seurat_label_f")], useNA="ifany")
#.........seurat_label_f
#cond_sex Micro.Vasc Astro L2.3 L4 Inhb L5 L6 Oligo
#...NTC F         18    20   20 20   20 20 20    17
#...NTC M         20    20   20 20   20 20 20    18
#...MDD F         20    20   20 20   20 20 20    16
#...MDD M         17    19   19 19   19 19 19    17
#...BPD F         20    20   20 20   20 20 20    18
#...BPD M         19    20   20 20   20 20 20    19


sil.results <- as.data.frame(approxSilhouette(reducedDim(spe_pseudo, "PCA_1663"), clusters=spe_pseudo$seurat_label_f))
sil.results$closest <- factor(ifelse(sil.results$width > 0, as.character(sil.results$cluster), as.character(sil.results$other)))
sil.results$closest <- factor(sil.results$closest, levels=levels(spe_pseudo$seurat_label_f))

p3 <- ggplot(sil.results, aes(x=cluster, y=width, colour=closest))+
  ggbeeswarm::geom_quasirandom()+scale_color_manual(values=cpList$transfer.bright)+
  labs(x="Seurat cell-type label", title="Seurat labels", subtitle="with low UMI cluster spots")+
  theme_minimal()+theme(text=element_text(size=10))

spe_se_low.umi = spe_pseudo

#label transfer (no low UMI)
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-no-lowUMI_sample-seurat-pc30_norm-filt.Rdata")
spe_pseudo$seurat_label_f = factor(as.character(spe_pseudo$seurat_label), 
                                   levels=c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"),
                                   labels=c("M.V","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))
table(spe_pseudo$seurat_label_f, useNA="ifany")
#Micro.Vasc      Astro       L2.3         L4       Inhb         L5         L6 
#114        119        119        119        119        119        119 
#Oligo 
#110 
spe_pseudo$cond_sex = factor(paste(spe_pseudo$condition, spe_pseudo$sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"))
table(as.data.frame(colData(spe_pseudo))[,c("cond_sex","seurat_label_f")], useNA="ifany")
#.........seurat_label_f
#cond_sex Micro.Vasc Astro L2.3 L4 Inhb L5 L6 Oligo
#...NTC F         18    20   20 20   20 20 20    18
#...NTC M         20    20   20 20   20 20 20    18
#...MDD F         20    20   20 20   20 20 20    17
#...MDD M         17    19   19 19   19 19 19    18
#...BPD F         20    20   20 20   20 20 20    19
#...BPD M         19    20   20 20   20 20 20    20

sil.results <- as.data.frame(approxSilhouette(reducedDim(spe_pseudo, "PCA_1663"), clusters=spe_pseudo$seurat_label_f))
sil.results$closest <- factor(ifelse(sil.results$width > 0, as.character(sil.results$cluster), as.character(sil.results$other)))
sil.results$closest <- factor(sil.results$closest, levels=levels(spe_pseudo$seurat_label_f))

p4 <- ggplot(sil.results, aes(x=cluster, y=width, colour=closest))+
  ggbeeswarm::geom_quasirandom()+scale_color_manual(values=cpList$transfer.bright)+
  labs(x="Seurat cell-type label", title="Seurat labels", subtitle="without low UMI cluster spots")+
  theme_minimal()+theme(text=element_text(size=10))

spe_se_no.low.umi = spe_pseudo


ggsave(file="plots/06_pseudobulk/low-UMI-cluster_silhouettes.png", 
       grid.arrange(p1, p2, p3, p4, ncol=2),
       bg="white", height=8, width=8)
cat("\nSaved silhouette plots to: plots/06_pseudobulk/low-UMI-cluster_silhouettes.png\n")

cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
