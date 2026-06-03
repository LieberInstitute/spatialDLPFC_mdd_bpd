setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(dplyr)
#  library(pheatmap)
  library(ggplot2)
  library(gridExtra)
#  library(ggrastr)
})

cpList <- readRDS("plots/colorPalettes.rds")
refined.modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv")

# aucell violins ----
aucell = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_modules-DEG-subset-refined_AUCell.csv", row.names=1)
colnames(aucell) = gsub("Regulon\\.for\\.","",colnames(aucell))

# normalize AUCell to help with upper limit
m1 = as.matrix(aucell[,1:(grep("seurat_label", colnames(aucell))-1)])

#scale aucell
## new conditional scale functions
find_3MAD = function(distribution) {
  median(distribution)+(3*mad(distribution))
}

find_q99 = function(distribution) {
  quantile(distribution, probs=.995)[[1]]
}

m1 = apply(m1, MARGIN=2, function(x) {
  thresh1 = find_3MAD(x)
  thresh2 = find_q99(x)
  nmax1 = sum(x>thresh1)
  nmax2 = sum(x>thresh2)
  if(nmax2<nmax1) {
    scale_max = thresh2
  } else {
    scale_max = thresh1
  }
  
  # if threshold>max, set threshold to max
  if(scale_max>max(x)) {
    scale_max = max(x)
  }
  
  # normalize to set max
  norm1 = x/scale_max
  norm1[norm1>1] = 1
  return(norm1)
  
})

aucell = cbind(as.data.frame(m1), aucell[,grep("seurat_label", colnames(aucell)):ncol(aucell)])

# now refined module subset (remove COX4I1 for redundancy and ADAMTS1 because too small)
mod_subset = c("A2M","IFITM3","CD74","HSPA1A","MT1X",
               "SNHG14","GLUL","CAMK2N1","GAD1","GRIN1","PRKAR1A",
               "UQCRH","APLP1","FTL","PLP1")

seurat_levels= c("M.V","Ast","L2.3","L4","Inb","L5","L6","Olg")
cond_sex = c("F NTC","F MDD", "F BPD","M NTC","M MDD","M BPD")
aucell$seurat_label= factor(aucell$seurat_label, levels=c("Micro.Vasc", "Astro", "L2.3", "L4", "Inhb", "L5", "L6", "Oligo"),
                            labels=seurat_levels)
aucell$sex.group = paste(aucell$sex, aucell$condition)
aucell$sex.group.se = paste(aucell$sex.group, aucell$seurat_label)
aucell$smoothed = factor(aucell$smoothed_k9_1663, levels=c("L1","L2","L3.4","L5","L6","WM","Vasc","GABA"),
                         labels=c("L1","L2","L3.4","L5","L6","WM","dropped","dropped"))
aucell$sex.group.sm = paste(aucell$sex.group, aucell$smoothed)

col.pal = cpList$transfer.bright[c(1:4,8,5:7)]
names(col.pal) = seurat_levels

aucell_long = tidyr::pivot_longer(aucell, all_of(mod_subset), names_to="module", values_to="AUCell") %>%
  mutate(module= factor(module, levels=mod_subset))


p1 <- ggplot(filter(aucell_long, module %in% c("PLP1","FTL","APLP1","HSPA1A")) %>%
	mutate(module=factor(module, levels=c("PLP1","FTL","APLP1","HSPA1A"))),
             aes(x=AUCell, color=seurat_label, lty=sex, group=sex.group.se))+
  stat_ecdf(linewidth=.5)+facet_grid(col=vars(module))+
  scale_x_continuous(breaks=c(0,.5,1), labels=c("0","0.5","1"))+
  scale_color_manual(values=col.pal, guide="none")+
  labs(title="Seurat label", x="AUCell (norm.)")+
  theme_minimal()+theme(axis.text.x=element_text(size=8),
	panel.grid.minor.y=element_blank())# , plot.margin= margin(6,16,6,6, "pt"))


p2 <- ggplot(filter(aucell_long, module %in% c("PLP1","FTL","APLP1","HSPA1A"), smoothed!="dropped") %>%
        mutate(module=factor(module, levels=c("PLP1","FTL","APLP1","HSPA1A"))),
             aes(x=AUCell, color=smoothed, lty=sex, group=sex.group.sm))+
  stat_ecdf(linewidth=.5)+facet_grid(col=vars(module))+
  scale_x_continuous(breaks=c(0,.5,1), labels=c("0","0.5","1"))+
  scale_color_manual(values=cpList$smoothed.bright, guide="none")+
  labs(title="PRECAST domain", x="AUCell (norm.)")+
  theme_minimal()+theme(axis.text.x=element_text(size=8),
	panel.grid.minor.y=element_blank())# , plot.margin= margin(6,16,6,6, "pt"))



p3 <- ggplot(filter(aucell_long, module %in% c("MT1X","IFITM3","A2M","CD74")) %>%
		mutate(module= factor(module, levels=c("MT1X","IFITM3","A2M","CD74"))),
             aes(x=AUCell, color=seurat_label, lty=sex, group=sex.group.se))+
  stat_ecdf(linewidth=.5)+facet_grid(col=vars(module))+
  scale_x_continuous(breaks=c(0,.5,1), labels=c("0","0.5","1"))+
  scale_color_manual(values=col.pal, guide="none")+
  labs(title="Seurat label", x="AUCell (norm.)")+
  theme_minimal()+theme(axis.text.x=element_text(size=8),
        panel.grid.minor.y=element_blank())# , plot.margin= margin(6,16,6,6, "pt"))


p4 <- ggplot(filter(aucell_long, module %in% c("MT1X","IFITM3","A2M","CD74"), smoothed!="dropped") %>%
		mutate(module=factor(module, levels=c("MT1X","IFITM3","A2M","CD74"))),
             aes(x=AUCell, color=smoothed, lty=sex, group=sex.group.sm))+
  stat_ecdf(linewidth=.5)+facet_grid(col=vars(module))+
  scale_x_continuous(breaks=c(0,.5,1), labels=c("0","0.5","1"))+
  scale_color_manual(values=cpList$smoothed.bright, guide="none")+
  labs(title="PRECAST domain", x="AUCell (norm.)")+
  theme_minimal()+theme(axis.text.x=element_text(size=8),
        panel.grid.minor.y=element_blank())# , plot.margin= margin(6,16,6,6, "pt"))

ggsave(file="plots/publication/Figure3/supp_select-modules-AUCell_ecdf.pdf",
	grid.arrange(p1, p2, p3, p4, ncol=1),
	width=7, height=8)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
