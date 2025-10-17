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

#Seurat PC20
cat("\nSeurat PC20 pseudobulk, before filtering...\n")
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc20_norm.Rdata")
spe_pseudo$seurat_label = factor(spe_pseudo$seurat_label, levels=c("Micro.Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"),
	labels=c("M.V","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"))
#spe_pseudo$seurat_label = factor(spe_pseudo$seurat_label, levels=c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"),
#        labels=c("M.V","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))
table(spe_pseudo$seurat_label, useNA="ifany")

spe_pseudo$cond_sex = factor(paste(spe_pseudo$condition, spe_pseudo$sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"))
table(as.data.frame(colData(spe_pseudo))[,c("cond_sex","seurat_label")], useNA="ifany")

sil.results <- as.data.frame(approxSilhouette(reducedDim(spe_pseudo, "PCA_1663"), clusters=spe_pseudo$seurat_label))
sil.results$closest <- factor(ifelse(sil.results$width > 0, as.character(sil.results$cluster), as.character(sil.results$other)))
sil.results$closest <- factor(sil.results$closest, levels=levels(spe_pseudo$seurat_label))

p1 <- ggplot(sil.results, aes(x=cluster, y=width, colour=closest))+
  ggbeeswarm::geom_quasirandom()+scale_color_manual(values=cpList$low.res.bright)+
  labs(x="Seurat label", title="Seurat PC20", subtitle="before QC filter")+
  theme_minimal()+theme(text=element_text(size=10))


#Seurat PC20 filtered
cat("\nSeurat PC20 pseudobulk, after filtering...\n")
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc20_norm-filt.Rdata")
spe_pseudo$seurat_label = factor(spe_pseudo$seurat_label, levels=c("Micro.Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"),
        labels=c("M.V","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"))
#spe_pseudo$seurat_label = factor(spe_pseudo$seurat_label, levels=c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"),
#        labels=c("M.V","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))
table(spe_pseudo$seurat_label, useNA="ifany")

spe_pseudo$cond_sex = factor(paste(spe_pseudo$condition, spe_pseudo$sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"))
table(as.data.frame(colData(spe_pseudo))[,c("cond_sex","seurat_label")], useNA="ifany")

sil.results <- as.data.frame(approxSilhouette(reducedDim(spe_pseudo, "PCA_1663"), clusters=spe_pseudo$seurat_label))
sil.results$closest <- factor(ifelse(sil.results$width > 0, as.character(sil.results$cluster), as.character(sil.results$other)))
sil.results$closest <- factor(sil.results$closest, levels=levels(spe_pseudo$seurat_label))

p2 <- ggplot(sil.results, aes(x=cluster, y=width, colour=closest))+
  ggbeeswarm::geom_quasirandom()+scale_color_manual(values=cpList$low.res.bright)+
  labs(x="Seurat label", title="Seurat PC20", subtitle="after QC filter")+
  theme_minimal()+theme(text=element_text(size=10))


#Seurat PC30
cat("\nSeurat PC30 pseudobulk, before filtering...\n")
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm.Rdata")
spe_pseudo$seurat_label	= factor(spe_pseudo$seurat_label, levels=c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"),
        labels=c("M.V","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))

table(spe_pseudo$seurat_label, useNA="ifany")

spe_pseudo$cond_sex = factor(paste(spe_pseudo$condition, spe_pseudo$sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"))
table(as.data.frame(colData(spe_pseudo))[,c("cond_sex","seurat_label")], useNA="ifany")


sil.results <- as.data.frame(approxSilhouette(reducedDim(spe_pseudo, "PCA_1663"), clusters=spe_pseudo$seurat_label))
sil.results$closest <- factor(ifelse(sil.results$width > 0, as.character(sil.results$cluster), as.character(sil.results$other)))
sil.results$closest <- factor(sil.results$closest, levels=levels(spe_pseudo$seurat_label))

p3 <- ggplot(sil.results, aes(x=cluster, y=width, colour=closest))+
  ggbeeswarm::geom_quasirandom()+scale_color_manual(values=cpList$transfer.bright)+
  labs(x="Seurat label", title="Seurat PC30", subtitle="before QC filter")+
  theme_minimal()+theme(text=element_text(size=10))


#Seurat PC30 filtered
cat("\nSeurat PC30 pseudobulk, after filtering...\n")
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
spe_pseudo$seurat_label = factor(spe_pseudo$seurat_label, levels=c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"),
        labels=c("M.V","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))

table(spe_pseudo$seurat_label, useNA="ifany")

spe_pseudo$cond_sex = factor(paste(spe_pseudo$condition, spe_pseudo$sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"))
table(as.data.frame(colData(spe_pseudo))[,c("cond_sex","seurat_label")], useNA="ifany")


sil.results <- as.data.frame(approxSilhouette(reducedDim(spe_pseudo, "PCA_1663"), clusters=spe_pseudo$seurat_label))
sil.results$closest <- factor(ifelse(sil.results$width > 0, as.character(sil.results$cluster), as.character(sil.results$other)))
sil.results$closest <- factor(sil.results$closest, levels=levels(spe_pseudo$seurat_label))

p4 <- ggplot(sil.results, aes(x=cluster, y=width, colour=closest))+
  ggbeeswarm::geom_quasirandom()+scale_color_manual(values=cpList$transfer.bright)+
  labs(x="Seurat label", title="Seurat PC30", subtitle="after QC filter")+
  theme_minimal()+theme(text=element_text(size=10))

# Seurat PC30 with no low UMI
cat("\nSeurat PC30 pseudobulked without low UMI spots...\n")
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-no-lowUMI_sample-seurat-pc30_norm.Rdata")
spe_pseudo$seurat_label = factor(spe_pseudo$seurat_label, levels=c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"),
        labels=c("M.V","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))
table(spe_pseudo$seurat_label, useNA="ifany")

spe_pseudo$cond_sex = factor(paste(spe_pseudo$condition, spe_pseudo$sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"))
table(as.data.frame(colData(spe_pseudo))[,c("cond_sex","seurat_label")], useNA="ifany")

sil.results <- as.data.frame(approxSilhouette(reducedDim(spe_pseudo, "PCA_1663"), clusters=spe_pseudo$seurat_label))
sil.results$closest <- factor(ifelse(sil.results$width > 0, as.character(sil.results$cluster), as.character(sil.results$other)))
sil.results$closest <- factor(sil.results$closest, levels=levels(spe_pseudo$seurat_label))

p5 <- ggplot(sil.results, aes(x=cluster, y=width, colour=closest))+
  ggbeeswarm::geom_quasirandom()+scale_color_manual(values=cpList$transfer.bright)+
  labs(x="Seurat label", title="Seurat PC30 (no low UMI spots)", subtitle="before QC filter")+
  theme_minimal()+theme(text=element_text(size=10))

ggsave(file="plots/06_pseudobulk/Seurat/Seurat_silhouettes.png", 
       grid.arrange(p1, p2, p3, p4, p3, p5, ncol=2),
       bg="white", height=10, width=8)
cat("\nSaved silhouette plots to: plots/06_pseudobulk/Seurat/Seurat_silhouettes.png\n")

cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
