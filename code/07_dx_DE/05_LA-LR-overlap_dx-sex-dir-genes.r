setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(dplyr)
	library(ggplot2)
	library(scater)
	library(pheatmap)
	library(gridExtra)
	library(grid)
	library(ggpubr)
})

set.seed(123)

cpList = readRDS("plots/colorPalettes.rds")

# load DE model results
## PRECAST smoothed
adj.results_sm = read.csv("processed-data/07_dx_DE/layer-adjusted-age_smoothed-k9-1663_compiled-results.csv", row.names = 1) %>%
  mutate(sex= factor(sex, levels=c("F","M")),
         group= factor(group, levels=c("NTC.MDD","NTC.BPD","MDD.BPD")),
         dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased", "increased")))

restr.results_sm <- read.csv("processed-data/07_dx_DE/layer-restricted-age_smoothed-k9-1663_compiled-results.csv", row.names=1) %>%
  mutate(sex= factor(sex, levels=c("F","M")),
         group= factor(group, levels=c("NTC.MDD","NTC.BPD","MDD.BPD")),
         smoothed=factor(cluster, levels=c("L1","L2","L3.4","L5","L6","WM")),
         dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased", "increased")))


## Seurat label transfer
adj.results_se = read.csv("processed-data/07_dx_DE/layer-adjusted-age_seurat-pc30-no-lowUMI_compiled-results.csv", row.names=1) %>%
  mutate(sex= factor(sex, levels=c("F","M")),
         group= factor(group, levels=c("NTC.MDD","NTC.BPD","MDD.BPD")),
         dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased", "increased")))

restr.results_se <- read.csv("processed-data/07_dx_DE/layer-restricted-age_seurat-pc30-no-lowUMI_compiled-results.csv", row.names=1) %>%
  mutate(sex= factor(sex, levels=c("F","M")),
         group= factor(group, levels=c("NTC.MDD","NTC.BPD","MDD.BPD")),
         seurat_label_f=factor(cluster, levels=c("Micro.Vasc","Astro","L2.3","L4","L5","L6","Oligo","Inhb"),
                               labels=c("M.V","Astro","L2.3","L4","L5","L6","Oligo","Inhb")),
         dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased", "increased")))

# load overlaps gene lists

saveList <- readRDS("processed-data/07_dx_DE/LA-LR-overlap_lists.rds")

# load spe for heatmap

load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
spe_sm <- spe_pseudo

load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-no-lowUMI_sample-seurat-pc30_norm-filt.Rdata")
spe_se <- spe_pseudo
spe_se$seurat_label_f = factor(spe_se$seurat_label, levels=c("Micro.Vasc","Astro","L2.3","L4","L5","L6","Inhb","Oligo"),
                               labels=c("M.V","Astro","L2.3","L4","L5","L6","Inhb","Oligo"))


# define dx*sex dir to plot

dx_sex_dir = "NTC.BPD_M_dn"


tmp = unlist(strsplit(dx_sex_dir, "_"))
target_group = tmp[[1]]
target_sex = tmp[[2]]
target_dir = tmp[[3]]
target_dir2 = ifelse(target_dir=="dn", "decreased", "increased")

# heatmaps
phm_test = plotGroupedHeatmap(spe_se, features = saveList[["LA.LR_both.annotations"]][[dx_sex_dir]], swap_rownames = "gene_name",
                              exprs_values="logcounts", group="seurat_label_f", 
                              center=T, fontsize_row=7,
                              angle_col=0, treeheight_row=10,
                              cluster_cols=F, clustering_method="ward.D2",
                              silent=T, 
                              main=paste0("Seurat: ", target_group," ", target_sex, ", ", target_dir2, "\n(logcount expr heatmap)"))

test_genes = phm_test$tree_row$labels[phm_test$tree_row$order]

phm_test2 = plotGroupedHeatmap(spe_sm, features= test_genes, swap_rownames="gene_name",
                               exprs_values="logcounts", group="smoothed_k9_1663",
                               center=T, fontsize_row=7,
                               angle_col=0,
                               cluster_cols=F, cluster_rows=F,
                               silent=T, 
                               main=paste0("PRECAST: ", target_group," ", target_sex, ", ", target_dir2, "\n(logcount expr heatmap)"))


plotList = list(arrangeGrob(phm_test2[[4]], phm_test[[4]], ncol=2))

# gene count table

t1 = as.data.frame(list("set"=c("L-A","L-A & L-R"),
                        "sm"=c(length(saveList[["sig_genes"]][[dx_sex_dir]]$adj_sm), length(saveList[["LA.LR_overlap"]][[dx_sex_dir]]$sm_both)),
                        "se"=c(length(saveList[["sig_genes"]][[dx_sex_dir]]$adj_se), length(saveList[["LA.LR_overlap"]][[dx_sex_dir]]$se_both)))
)
gt1 = tableGrob(t1, rows = NULL, cols=c(paste("Sig.", target_dir2, sep="\n"),
                                        "PRECAST\n(smoothed)","Seurat labels"))

if(target_dir=="dn") fill_color = "skyblue"
if(target_dir=="up") fill_color = "tomato"
t2 = as.data.frame(list("set"="L-A & L-R\n(both annotations)",
                        "sm"=length(saveList[["LA.LR_both.annotations"]][[dx_sex_dir]])))
tt2 <- ttheme_default(core=list(bg_params = list(fill=c(fill_color,fill_color)))
)
gt2 = tableGrob(t2, rows = NULL, cols=NULL, theme=tt2)

gt3 = gtable_combine(gt1, gt2, along=2, join="left")
gt3$layout[c(20,22),"r"] = 3


# L-R p sig legend

gt = as.data.frame(list(symbol=c("","^","*","**","***"),
                        adj.P.Val=c("> 0.1","< 0.1","< 0.05","< 0.01","< 0.001"))
)
plotList[[2]] = arrangeGrob(gt3, tableGrob(gt, rows = NULL), 
                            layout_matrix=t(matrix(c(NA,1,1,1,NA,2,2,NA))), 
                            top="Left: overlap breakdown; Right: sig. key")



# collect and format DF for barplots

## precast smoothed
spe_tmp = spe_sm[,spe_sm$sex==target_sex & spe_sm$condition %in% unlist(strsplit(target_group, "\\."))]

for (j in test_genes) {
  colData(spe_tmp)[[j]] = logcounts(spe_tmp)[rowData(spe_tmp)$gene_name==j,]
}

#preserve names so if any names have '-' character they don't throw an error
preserve.names = gsub("-", "\\.", test_genes)

plot.df = as.data.frame(colData(spe_tmp)[,c("condition", "sex", "smoothed_k9_1663", "sample_id", test_genes)]) %>%
  tidyr::pivot_longer(all_of(preserve.names), names_to="key_genes", values_to="logcounts") %>%
  mutate(key_genes= factor(key_genes, levels=preserve.names, labels=test_genes),
         x_labels=as.character(smoothed_k9_1663))


## seurat labels
spe_tmp2 = spe_se[,spe_se$sex==target_sex & spe_se$condition %in% unlist(strsplit(target_group, "\\."))]

for (j in test_genes) {
  colData(spe_tmp2)[[j]] = logcounts(spe_tmp2)[rowData(spe_tmp2)$gene_name==j,]
}

plot.df2 = as.data.frame(colData(spe_tmp2)[,c("condition", "sex", "seurat_label_f", "sample_id", test_genes)]) %>%
  tidyr::pivot_longer(all_of(preserve.names), names_to="key_genes", values_to="logcounts") %>%
  mutate(key_genes= factor(key_genes, levels=preserve.names, labels=test_genes),
         x_labels= as.character(seurat_label_f))

## combine

plot.both = bind_rows(mutate(plot.df[,intersect(colnames(plot.df), colnames(plot.df2))], annotation="sm"),
                      mutate(plot.df2[,intersect(colnames(plot.df), colnames(plot.df2))], annotation="se")) %>%
  mutate(annotation = factor(annotation, levels=c("sm","se"), labels=c("PRECAST (smoothed)", "Seurat labels")),
         x_labels= factor(x_labels, levels=c("M.V","Astro","L1","L2.3","L2","L3.4","L4","L5","L6","Inhb","WM","Oligo")))


# generate barplot list

plist <- lapply(test_genes, function(x) {
  suppressMessages({
    tmp = filter(plot.both, key_genes==x)
    ylim1 = c(min(tmp$logcounts), max(tmp$logcounts))
    ylim2 = ylim1 #copy for editing and comparing
    
    la.stat = filter(adj.results_sm, group==target_group, sex==target_sex, gene_name==x)
    subt = paste("L-A sig: adj p<", format(la.stat$adj.P.Val, scientific=T, digits=1),
                 "logFC =", round(la.stat$logFC, 2))
    bp1 = ggplot(filter(tmp, annotation=="PRECAST (smoothed)"), 
                 aes(x=x_labels, y=logcounts, fill=condition))+
      geom_boxplot(position=position_dodge(width=.6), width=.5, outliers = F)+
      scale_fill_manual(values=cpList$dx.pal)+
      labs(title="PRECAST (smoothed)", subtitle=subt)+ylim(ylim1)+
      theme_bw()+theme(legend.position="bottom", legend.title = element_blank(), axis.title.x=element_blank())
    gbp1 = ggplot_build(bp1)$data[[1]]
    gbp1$x_labels = rep(levels(droplevels(filter(tmp, annotation=="PRECAST (smoothed)")$x_labels)), each=2)

    tmp1 = filter(restr.results_sm, gene_name==x, group==target_group, sex==target_sex) %>%
      mutate(group1=unlist(strsplit(target_group, "\\."))[[1]], group2=unlist(strsplit(target_group, "\\."))[[2]],
             p.signif = cut(adj.P.Val, breaks=c(0, .001, .01, .05, .1), labels=c("***","**","*","^"))
      )
    tmp.stats = filter(plot.both, key_genes==x, annotation=="PRECAST (smoothed)", condition %in% unlist(strsplit(target_group, "\\.")), sex==target_sex) %>%
      group_by(condition, sex, x_labels) %>% 
      summarise(n=n(), med_test=median(logcounts),
	nup= med_test+1.58*IQR(logcounts)/sqrt(n)) %>%
      left_join(gbp1, by=c("x_labels","med_test"="middle","nup"="notchupper")) %>%
      group_by(x_labels) %>% summarise(max_ymax= max(ymax), more_ymax= max_ymax*1.01)
    stat.test = left_join(tmp1[,c("group1","group2","gene_name","sex","smoothed","adj.P.Val","p.signif")],
                          tmp.stats, by=c("smoothed"="x_labels"))
    if(max(stat.test$more_ymax) > ylim2[[2]]) {
      ylim2[[2]] = max(stat.test$more_ymax) 
    }
    bp1 = bp1 + stat_pvalue_manual(filter(stat.test, !is.na(p.signif)), x="smoothed", 
                                   y.position="more_ymax", 
                                   label="p.signif", position=position_dodge(width=.6),
                                   size=6)
    
    la.stat = filter(adj.results_se, group==target_group, sex==target_sex, gene_name==x)
    subt = paste("L-A sig: adj p<", format(la.stat$adj.P.Val, scientific=T, digits=1),
                 "logFC =", round(la.stat$logFC, 2))
    bp2 = ggplot(filter(tmp, annotation=="Seurat labels"), 
                 aes(x=x_labels, y=logcounts, fill=condition))+
      geom_boxplot(position=position_dodge(width=.6), width=.5, outliers = F)+
      scale_fill_manual(values=cpList$dx.pal)+
      labs(title="Seurat labels", subtitle=subt)+ylim(ylim1)+
      theme_bw()+theme(legend.position="bottom", legend.title = element_blank(), axis.title.x=element_blank())
    gbp2 = ggplot_build(bp2)$data[[1]]
    gbp2$x_labels = rep(levels(droplevels(filter(tmp, annotation=="Seurat labels")$x_labels)), each=2)

    tmp2 = filter(restr.results_se, gene_name==x, group==target_group, sex==target_sex) %>%
      mutate(group1=unlist(strsplit(target_group, "\\."))[[1]], group2=unlist(strsplit(target_group, "\\."))[[2]],
             p.signif = cut(adj.P.Val, breaks=c(0, .001, .01, .05, .1), labels=c("***","**","*","^")),
             seurat_label_f= factor(seurat_label_f, levels=levels(tmp$x_labels))
      )
    tmp.stats = filter(plot.both, key_genes==x, annotation=="Seurat labels", condition %in% unlist(strsplit(target_group, "\\.")), sex==target_sex) %>%
      group_by(condition, sex, x_labels) %>% 
      summarise(n=n(), med_test=median(logcounts),
	nup= med_test+1.58*IQR(logcounts)/sqrt(n)) %>%
      left_join(gbp2, by=c("x_labels","med_test"="middle","nup"="notchupper")) %>%
      group_by(x_labels) %>% summarise(max_ymax= max(ymax), more_ymax= max_ymax*1.01)
    stat.test = left_join(tmp2[,c("group1","group2","gene_name","sex","seurat_label_f","adj.P.Val","p.signif")],
                          tmp.stats, by=c("seurat_label_f"="x_labels"))
    if(max(stat.test$more_ymax) > ylim2[[2]]) {
      ylim2[[2]] = max(stat.test$more_ymax) 
    }
    bp2 = bp2 + stat_pvalue_manual(filter(stat.test, !is.na(p.signif)), x="seurat_label_f", 
                                   y.position="more_ymax", 
                                   label="p.signif", position=position_dodge(width=.6),
                                   size=6)
    
  })
  
  if(identical(ylim1, ylim2)) {
    arrangeGrob(bp1, bp2, ncol=2, top=textGrob(x, gp=gpar(fontsize=14,fontface=4)))
  } else {
    arrangeGrob(bp1+ylim(ylim2), bp2+ylim(ylim2), ncol=2, top=textGrob(x, gp=gpar(fontsize=14,fontface=4)))
  }
  
})

##for NTC.MDD F up, because there are so many, make 2 rows per page and make pages bigger
#plist = marrangeGrob(plist, nrow=2, ncol=1, top = NULL)

plotList = c(plotList, plist)


# save plots to multi-page pdf

ggsave(file=paste0("plots/07_dx_DE/LA-LR-overlap_", target_group, "-", target_sex, "-", target_dir2, ".pdf"),
       marrangeGrob(plotList, nrow=1, ncol=1, top = NULL),
       height=5, width=8)
	#height=11, width=8)
cat("\n\nSaved compiled pdf to:", paste0("plots/07_dx_DE/LA-LR-overlap_", target_group, "-", target_sex, "-", target_dir2, ".pdf"))


cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
