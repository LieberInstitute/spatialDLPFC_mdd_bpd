setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
})

aucell = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_modules-DEG-subset-refined_AUCell.csv")
colnames(aucell) = gsub("Regulon\\.for\\.","",colnames(aucell))
aucell$sample_id = substr(rownames(aucell), start=20, stop=50)
aucell$seurat_label = factor(aucell$seurat_label, levels=c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"),
	labels=c("M/V","Ast","L2/3","L4","Inb","L5","L6","Olg"))

tmp <- group_by(aucell, seurat_label, sample_id) %>% summarise(n=n(), gad1.pos=sum(GAD1>0), prop.gad1.pos=gad1.pos/n, avg.gad1=mean(GAD1), med.gad1=median(GAD1))
tmp$avg.gad1_scale = scale(tmp$avg.gad1)

p1 <- ggplot(tmp, aes(x=seurat_label, y=prop.gad1.pos, color=avg.gad1_scale))+
  ggbeeswarm::geom_quasirandom(size=.5)+
  scale_color_gradient("Donor avg.\nGAD1-module\nAUCell (z-score)",low="lightgrey", high="black")+
  scale_y_continuous(limits=c(0,1), breaks=c(0,.5,1))+labs(x="domain-CT", y="prop. spots GAD1-module > 0")+
  theme_minimal()+theme(text=element_text(size=6), legend.key.width=unit(8,"pt"), legend.key.height=unit(10,"pt"))

ggsave(file="plots/publication/Figure_nrn/GAD1-module_by-domain-CT.pdf", p1, height=3, width=3.5)

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
