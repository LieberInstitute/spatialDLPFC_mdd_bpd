setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(pheatmap)
})

set.seed(123)

mod_subset = c("SNHG14","PRKAR1A","COX4I1","EEF1A1",
  "PLP1",
  "GFAP","GLUL","HSPA1A","MT1M",
  "IFITM3","CD74","A2M",
  "GRIN1","CAMK2N1","GAD1")

aucell = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_modules-DEG-subset-refined_AUCell.csv")
colnames(aucell) = gsub("Regulon\\.for\\.","",colnames(aucell))

seurat_levels= c("M.V","Ast","L2.3","L4","Inb","L5","L6","Olg")
cond_sex = c("F NTC","F MDD", "F BPD","M NTC","M MDD","M BPD")

aucell$seurat_label= factor(aucell$seurat_label, levels=c("Micro.Vasc", "Astro", "L2.3", "L4", "Inhb", "L5", "L6", "Oligo"),
	labels=seurat_levels)
aucell$sex.group = paste(aucell$sex, aucell$condition)
aucell$clus_groups= factor(paste(aucell$sex.group, aucell$seurat_label),
	levels=as.character(outer(cond_sex, seurat_levels, paste)))

# calculate averages per dx*sex*cluster
avg.df = group_by(aucell, clus_groups) %>% summarise_at(mod_subset, mean)
avg.mtx = as.matrix(avg.df[,-1])
rownames(avg.mtx) = avg.df$clus_groups

# annotate columns by cluster
col_annot = data.frame("cluster"=unlist(lapply(strsplit(as.character(avg.df$clus_groups), split = " "), function(x) x[[3]])), row.names = avg.df$clus_groups)
cpList <- readRDS("plots/colorPalettes.rds")
annot_colors = list("cluster"=cpList$transfer.bright[c(1:4,8,5:7)])
names(annot_colors$cluster) = seurat_levels

phm3 = pheatmap(t(avg.mtx), cluster_rows=T, cluster_cols=F, border_color=NA,
	clustering_method="ward.D2", treeheight_row=20,
         color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name ="Greys"))(100),
         annotation_col=col_annot, annotation_colors=annot_colors,
         annotation_legend=F, annotation_names_row=F, annotation_names_col=F,
         angle_col=90, fontsize_col = 6, scale="row",
         breaks= seq(from = -1, to = 3, length.out = 101))

# and now plot by smoothed domain for supplement
precast_levels= c("L1","L2","L3.4","L5","L6","WM")

avg.df2 = filter(aucell, smoothed_k9_1663 %in% precast_levels) %>% 
  mutate(smoothed_k9_1663=factor(smoothed_k9_1663, levels=precast_levels),
                clus_groups= factor(paste(sex.group, smoothed_k9_1663),
                                    levels=as.character(outer(cond_sex, precast_levels, paste)))) %>%
  group_by(clus_groups) %>% summarise_at(mod_subset, mean)

avg.mtx2 = as.matrix(avg.df2[,-1])
rownames(avg.mtx2) = avg.df2$clus_groups
col_annot2 = data.frame("cluster"=unlist(lapply(strsplit(as.character(avg.df2$clus_groups), split = " "), function(x) x[[3]])), row.names = avg.df2$clus_groups)
annot_colors = list("cluster"=cpList$smoothed.bright)

phm4 = pheatmap(t(avg.mtx2), cluster_rows=T, cluster_cols=F, border_color=NA,
	clustering_method="ward.D2", treeheight_row=20,
         color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name ="Greys"))(100),
         annotation_col=col_annot2, annotation_colors=annot_colors,
         annotation_legend=F, annotation_names_row=F, annotation_names_col=F,
         angle_col=90, fontsize_col = 6, scale="row",
         breaks= seq(from = -1, to = 3, length.out = 101))


pdf(file="plots/publication/modules/modules_avg-AUCell_heatmaps.pdf", width=5, height=5)
plot(phm3[[4]])
plot(phm4[[4]])
dev.off()

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
