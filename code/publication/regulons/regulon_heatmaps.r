setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(dplyr)
	library(pheatmap)
})
set.seed(123)


# no regulons smaller than 10
regulons <- read.csv("processed-data/10_SCENIC/spe-n109_13162-no-lowUMI_logcounts_regulons-top20_merged.csv", row.names=1)
reg_subset = filter(regulons, set_size>=10)$TF

aucell_corr = read.csv("processed-data/10_SCENIC/spe-n109_13162-no-lowUMI_regulons-top20_AUCell-all-spots-correlation.csv", row.names=1)
colnames(aucell_corr) = substr(colnames(aucell_corr), start=0, stop=nchar(colnames(aucell_corr))-3)

col.pal = rev(RColorBrewer::brewer.pal("RdBu", n=8))
col.pal = c(col.pal[1:4], "white", col.pal[5:8])

phm1 = pheatmap(aucell_corr[paste0(reg_subset, "(+)"),reg_subset],
        color = colorRampPalette(col.pal)(100), breaks = seq(-1, 1, length.out=101),
        clustering.method="ward.D2", angle_col=90,
        main="AUCell correlations")

# now actual AUCell values
aucell = read.csv("processed-data/10_SCENIC/spe-n109_13162-no-lowUMI_regulons-top20_AUCell.csv", row.names=1)
replace_names = 1:(grep("seurat_label", colnames(aucell))-1)
colnames(aucell)[replace_names] = substr(colnames(aucell)[replace_names], start=0, stop=nchar(colnames(aucell)[replace_names])-3)

## correlate with module aucell values
aucell_modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_modules-DEG-subset-refined_AUCell.csv", row.names=1)
colnames(aucell_modules) = gsub("Regulon\\.for\\.", "mod_", colnames(aucell_modules))
aucell_modules2 = aucell_modules[rownames(aucell),]

cor.mtx = cor(aucell_modules2[,grep("mod_", colnames(aucell_modules))], aucell[,reg_subset])
#mod_subset = c("A2M","IFITM3","CD74","HSPA1A","MT1X","SNHG14","GLUL",
#               "CAMK2N1","GAD1","GRIN1","PRKAR1A","UQCRH","APLP1","FTL","PLP1")
mod_subset = c("SNHG14","PRKAR1A","COX4I1","EEF1A1",
  "PLP1",
  "GFAP","GLUL","HSPA1A","MT1M",
  "IFITM3","CD74","A2M",
  "GRIN1","CAMK2N1","GAD1")

rownames(cor.mtx) = gsub("mod_","", rownames(cor.mtx))

## now make secondary matrix describing the overlap of genes to label with values
regList = strsplit(regulons$set_str, "/")
names(regList) = regulons$TF
#regList = regList[reg_subset]

refined.modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv")
modList = mod_subset
names(modList) = mod_subset
modList = lapply(modList, function(x) c(x, filter(refined.modules, TF==x)$target))

df1 = do.call(rbind, lapply(reg_subset, function(i) {
  data.frame("regulon"=i, "reg_size"=length(regList[[i]]), "module"=mod_subset, "mod_size"=sapply(modList, length),
             "intersect_size"=sapply(modList, function(x) length(intersect(x, regList[[i]]))),
             "jcoef"=sapply(modList, function(x) length(intersect(x, regList[[i]])))/sapply(modList, function(x) length(union(x, regList[[i]]))))
}))
df2 = tidyr::pivot_wider(df1[c("regulon","module","intersect_size")], names_from="module", values_from="intersect_size")
m2 = as.matrix(df2[,-1])
rownames(m2) = df2$regulon
m3 = m2
m3[m3==0] = ""

phm2 = pheatmap(t(cor.mtx[mod_subset,reg_subset]), #cluster_row=F, 
	clustering_method="ward.D2", treeheight_row=20, treeheight_col=20,
         color = colorRampPalette(col.pal)(100), breaks = seq(-1, 1, length.out=101),
	 display_numbers=m3[reg_subset,mod_subset], number_col="black",
         angle_col=90, main= "AUCell correlation (DEG modules vs regulons)")

# now avg heatmaps
seurat_levels= c("M.V","Ast","L2.3","L4","Inb","L5","L6","Olg")
cond_sex = c("F NTC","F MDD", "F BPD","M NTC","M MDD","M BPD")

aucell$seurat_label= factor(aucell$seurat_label, levels=c("Micro.Vasc", "Astro", "L2.3", "L4", "Inhb", "L5", "L6", "Oligo"),
	labels=seurat_levels)
aucell$sex.group = paste(aucell$sex, aucell$condition)
aucell$clus_groups= factor(paste(aucell$sex.group, aucell$seurat_label),
	levels=as.character(outer(cond_sex, seurat_levels, paste)))

# calculate averages per dx*sex*cluster
avg.df = group_by(aucell, clus_groups) %>% summarise_at(reg_subset, mean)
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
  group_by(clus_groups) %>% summarise_at(reg_subset, mean)

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


pdf(file="plots/publication/regulons/regulon_heatmaps.pdf", width=5, height=5)
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
