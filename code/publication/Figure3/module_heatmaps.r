setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(pheatmap)
})

set.seed(123)

# first plot all DEG modules
input.mtx = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_DEG-modules_target-pearson-correlation.csv", row.names=1)
# make dimnames consistent ("-" to ".")
rownames(input.mtx) = gsub("-", "\\.", rownames(input.mtx))

# top predictors (unfiltered)
mod_subset = c("COX4I1","UQCRH","PRKAR1A","CAMK2N1","GRIN1","GAD1",
               "SNHG14","GLUL","A2M","ADAMTS1","IFITM3","CD74",
               "FTL","APLP1","PLP1","HSPA1A","MT1X")

col_annot = data.frame("is_mod"=as.character(colnames(input.mtx) %in% mod_subset), row.names=colnames(input.mtx))
annot_colors = list("is_mod"=c("FALSE"="white", "TRUE"="black"))

phm1 = pheatmap(input.mtx, clustering_method="ward.D2", border_color=NA,
	breaks= seq(from = -1, to = 1, length.out = 101),
	treeheight_row=20, treeheight_col=20, show_rownames=F, show_colnames=F, 
	annotation_legend=F, annotation_names_row=F, annotation_names_col=F,
	annotation_row=col_annot, annotation_col=col_annot, annotation_colors=annot_colors)

# now plot top predictor DEG modules
phm2 = pheatmap(input.mtx[mod_subset, mod_subset], clustering_method="ward.D2",	border_color=NA,
        breaks= seq(from = -1, to = 1, length.out = 101),
        treeheight_row=20, treeheight_col=20, angle_col=90)

# now refined module subset (remove COX4I1 for redundancy and ADAMTS1 because too small)
mod_subset = c("A2M","IFITM3","CD74","HSPA1A","MT1X",
	"SNHG14","GLUL","CAMK2N1","GAD1","GRIN1","PRKAR1A",
	"UQCRH","APLP1","FTL","PLP1")

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

phm3 = pheatmap(t(avg.mtx[,mod_subset]), cluster_rows=F, cluster_cols=F, border_color=NA,
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

phm4 = pheatmap(t(avg.mtx2[,mod_subset]), cluster_rows=F, cluster_cols=F, border_color=NA,
         color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name ="Greys"))(100),
         annotation_col=col_annot2, annotation_colors=annot_colors,
         annotation_legend=F, annotation_names_row=F, annotation_names_col=F,
         angle_col=90, fontsize_col = 6, scale="row",
         breaks= seq(from = -1, to = 3, length.out = 101))


pdf(file="plots/publication/Figure3/module_heatmaps.pdf", width=5, height=5)
plot(phm1[[4]])
plot(phm2[[4]])
plot(phm3[[4]])
plot(phm4[[4]])
dev.off()

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()

