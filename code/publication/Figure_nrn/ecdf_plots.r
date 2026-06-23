setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
})

cpList <- readRDS("plots/colorPalettes.rds")

# from modules
#aucell = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_modules-DEG-subset-refined_AUCell.csv", row.names=1)
#colnames(aucell) = gsub("Regulon\\.for\\.","",colnames(aucell))
#seurat_levels= c("M.V","Ast","L2.3","L4","Inb","L5","L6","Olg")
#aucell$seurat_label= factor(aucell$seurat_label, levels=c("Micro.Vasc", "Astro", "L2.3", "L4", "Inhb", "L5", "L6", "Oligo"),
#                            labels=seurat_levels)

#col.pal = cpList$transfer.bright[c(1:4,8,5:7)]
#names(col.pal) = seurat_levels

# from regulons
## no regulons smaller than 10
regulons <- read.csv("processed-data/10_SCENIC/spe-n109_13162-no-lowUMI_logcounts_regulons-top20_merged.csv", row.names=1)
reg_subset = filter(regulons, set_size>=10)$TF
aucell = read.csv("processed-data/10_SCENIC/spe-n109_13162-no-lowUMI_regulons-top20_AUCell.csv", row.names=1)
replace_names = 1:(grep("seurat_label", colnames(aucell))-1)
colnames(aucell)[replace_names] = substr(colnames(aucell)[replace_names], start=0, stop=nchar(colnames(aucell)[replace_names])-3)
precast_levels= c("L1","L2","L3.4","L5","L6","WM")
aucell$smoothed_k9_1663 = factor(aucell$smoothed_k9_1663, levels=precast_levels)

col.pal = cpList$smoothed.bright

#mod_subset = c("GRIN1","CAMK2N1","GAD1")
mod_subset = c("ARX","LHX6","DLX1")

aucell_long = tidyr::pivot_longer(aucell, all_of(mod_subset), names_to="module", values_to="AUCell") %>%
  mutate(module= factor(module, levels=mod_subset))

p1 <- #ggplot(aucell_long, aes(x=AUCell, color=seurat_label))+#, lty=sex, group=sex.group.se))+
	ggplot(aucell_long, aes(x=AUCell, color=smoothed_k9_1663))+
  stat_ecdf()+facet_grid(col=vars(module), scales="free_x")+
  scale_y_continuous(breaks=c(0,.5,1), labels=c("0","0.5","1"))+
  scale_color_manual(values=col.pal, guide="none")+
  labs(x="AUCell")+
  theme_minimal()+theme(text=element_text(size=6),
                        panel.grid.minor.x=element_blank())

#ggsave(file="plots/publication/Figure_nrn/GRIN1-CAMK2N1-GAD1_ecdf.pdf", p1, height=1.7, width=3.75)
ggsave(file="plots/publication/Figure_nrn/ARX-LHX6-DLX1_ecdf.pdf", p1, height=1.7, width=3.75)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()

