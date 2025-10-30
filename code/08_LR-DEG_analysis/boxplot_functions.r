suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
	library(grid)
	library(ggpubr)
})

### load required data
loadData <- function() {
  load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
  spe_sm <- spe_pseudo
  
  load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
  spe_se <- spe_pseudo
  spe_se$seurat_label_f = factor(spe_se$seurat_label, levels=c("Micro.Vasc","Astro","L2.3","L4","L5","L6","Inhb","Oligo"),
                                 labels=c("M.V","Astro","L2.3","L4","L5","L6","Inhb","Oligo"))
  
  # collect and format DF for barplots
  # load DE model results
  ## PRECAST smoothed
  adj.results_sm = read.csv("processed-data/07_dx_DE/layer-adjusted-age_smoothed-k9-1663_rev-gene-input_compiled-results.csv", row.names = 1) %>%
    mutate(sex= factor(sex, levels=c("F","M")),
           group= factor(group, levels=c("NTC.MDD","NTC.BPD","MDD.BPD")),
           dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased", "increased")))
  
  restr.results_sm <- read.csv("processed-data/07_dx_DE/layer-restricted-age_smoothed-k9-1663_rev-gene-input_compiled-results.csv", row.names=1) %>%
    mutate(sex= factor(sex, levels=c("F","M")),
           group= factor(group, levels=c("NTC.MDD","NTC.BPD","MDD.BPD")),
           smoothed=factor(cluster, levels=c("L1","L2","L3.4","L5","L6","WM")),
           dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased", "increased")))
  
  
  ## Seurat label transfer
  adj.results_se = read.csv("processed-data/07_dx_DE/layer-adjusted-age_seurat-pc30_rev-gene-input_compiled-results.csv", row.names=1) %>%
    mutate(sex= factor(sex, levels=c("F","M")),
           group= factor(group, levels=c("NTC.MDD","NTC.BPD","MDD.BPD")),
           dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased", "increased")))
  
  restr.results_se <- read.csv("processed-data/07_dx_DE/layer-restricted-age_seurat-pc30_rev-gene-input_compiled-results.csv", row.names=1) %>%
    mutate(sex= factor(sex, levels=c("F","M")),
           group= factor(group, levels=c("NTC.MDD","NTC.BPD","MDD.BPD")),
           seurat_label_f=factor(cluster, levels=c("Micro.Vasc","Astro","L2.3","L4","L5","L6","Oligo","Inhb"),
                                 labels=c("M.V","Astro","L2.3","L4","L5","L6","Oligo","Inhb")),
           dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased", "increased")))
  
  return(list("spe_smoothed"=spe_sm, "LA_smoothed"=adj.results_sm, "LR_smoothed"=restr.results_sm,
              "spe_seurat"=spe_se, "LA_seurat"=adj.results_se, "LR_seurat"=restr.results_se
              ))
}

### take gene set and subsetted spe objects and generate DF formatted for boxplotting
formatData <- function(gene_set, spe_sm, spe_se) {
  ## precast smoothed
  for (j in gene_set) {
    colData(spe_sm)[[j]] = logcounts(spe_sm)[rowData(spe_sm)$gene_name==j,]
  }
  
  #preserve names so if any names have '-' character they don't throw an error
  preserve.names = gsub("-", "\\.", gene_set)
  
  plot.df = as.data.frame(colData(spe_sm)[,c("condition", "sex", "smoothed_k9_1663", "sample_id", gene_set)]) %>%
    tidyr::pivot_longer(all_of(preserve.names), names_to="key_genes", values_to="logcounts") %>%
    mutate(key_genes= factor(key_genes, levels=preserve.names, labels=gene_set),
           x_labels=as.character(smoothed_k9_1663))
  
  
  ## seurat labels
  for (j in gene_set) {
    colData(spe_se)[[j]] = logcounts(spe_se)[rowData(spe_se)$gene_name==j,]
  }
  
  plot.df2 = as.data.frame(colData(spe_se)[,c("condition", "sex", "seurat_label_f", "sample_id", gene_set)]) %>%
    tidyr::pivot_longer(all_of(preserve.names), names_to="key_genes", values_to="logcounts") %>%
    mutate(key_genes= factor(key_genes, levels=preserve.names, labels=gene_set),
           x_labels= as.character(seurat_label_f))
  
  
  plot.both = bind_rows(mutate(plot.df[,intersect(colnames(plot.df), colnames(plot.df2))], annotation="sm"),
                        mutate(plot.df2[,intersect(colnames(plot.df), colnames(plot.df2))], annotation="se")) %>%
    mutate(annotation = factor(annotation, levels=c("sm","se"), labels=c("PRECAST (smoothed)", "Seurat labels")),
           x_labels= factor(x_labels, levels=c("M.V","Astro","L1","L2.3","L2","L3.4","L4","Inhb","L5","L6","WM","Oligo")))
  
  return(list("spe_smoothed"=spe_sm, "spe_seurat"=spe_se, "boxplot.df"=plot.both))
}

### generate paired boxplot for both annotations
plotBoxplot <- function(plot_gene, plot_both, DE_list, plot_outliers=F) {
  adj.results_sm = DE_list[["LA_smoothed"]]
  adj.results_se = DE_list[["LA_seurat"]]
  restr.results_sm = DE_list[["LR_smoothed"]]
  restr.results_se = DE_list[["LR_seurat"]]

  suppressMessages({
    tmp = filter(plot_both, key_genes==plot_gene)
    ylim1 = c(min(tmp$logcounts), max(tmp$logcounts))
    ylim2 = ylim1 #copy for editing and comparing
    
    la.stat = filter(adj.results_sm, group==target_group, sex==target_sex, gene_name==plot_gene)
    if(nrow(la.stat)==0) {
        subt = "Gene not present in annotation pseudobulk"
        subt_color="grey50"
    } else {
        subt = paste("L-A sig: adj p<", format(la.stat$adj.P.Val, scientific=T, digits=1),
                 "logFC =", round(la.stat$logFC, 2))
        subt_color=ifelse(la.stat$adj.P.Val<.05, "black", "grey50")
    }
    bp1 = ggplot(filter(tmp, annotation=="PRECAST (smoothed)"), 
                 aes(x=x_labels, y=logcounts, fill=condition))+
      geom_boxplot(position=position_dodge(width=.6), width=.5, outliers = plot_outliers)+
      scale_fill_manual(values=cpList$dx.pal)+
      labs(title="PRECAST (smoothed)", subtitle=subt)+ylim(ylim1)+
      theme_bw()+theme(legend.position="bottom", legend.title = element_blank(), axis.title.x=element_blank(),
		plot.subtitle=element_text(color=subt_color))
    gbp1 = ggplot_build(bp1)$data[[1]]
    gbp1$x_labels = rep(levels(droplevels(filter(tmp, annotation=="PRECAST (smoothed)")$x_labels)), each=2)
    
    tmp1 = filter(restr.results_sm, gene_name==plot_gene, group==target_group, sex==target_sex) %>%
      mutate(group1=unlist(strsplit(target_group, "\\."))[[1]], group2=unlist(strsplit(target_group, "\\."))[[2]],
             p.signif = cut(adj.P.Val, breaks=c(0, .001, .01, .05, .1), labels=c("***","**","*","^"))
      )
    tmp.stats = filter(plot_both, key_genes==plot_gene, annotation=="PRECAST (smoothed)", condition %in% unlist(strsplit(target_group, "\\.")), sex==target_sex) %>%
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
    
    la.stat = filter(adj.results_se, group==target_group, sex==target_sex, gene_name==plot_gene)
    if(nrow(la.stat)==0) {
        subt = "Gene not present in annotation pseudobulk"
	subt_color="grey50"
    } else {
        subt = paste("L-A sig: adj p<", format(la.stat$adj.P.Val, scientific=T, digits=1),
                 "logFC =", round(la.stat$logFC, 2))
	subt_color=ifelse(la.stat$adj.P.Val<.05, "black", "grey50")
    }
    bp2 = ggplot(filter(tmp, annotation=="Seurat labels"), 
                 aes(x=x_labels, y=logcounts, fill=condition))+
      geom_boxplot(position=position_dodge(width=.6), width=.5, outliers = plot_outliers)+
      scale_fill_manual(values=cpList$dx.pal)+
      labs(title="Seurat labels", subtitle=subt)+ylim(ylim1)+
      theme_bw()+theme(legend.position="bottom", legend.title = element_blank(), axis.title.x=element_blank(),
		plot.subtitle=element_text(color=subt_color))
    gbp2 = ggplot_build(bp2)$data[[1]]
    gbp2$x_labels = rep(levels(droplevels(filter(tmp, annotation=="Seurat labels")$x_labels)), each=2)
    
    tmp2 = filter(restr.results_se, gene_name==plot_gene, group==target_group, sex==target_sex) %>%
      mutate(group1=unlist(strsplit(target_group, "\\."))[[1]], group2=unlist(strsplit(target_group, "\\."))[[2]],
             p.signif = cut(adj.P.Val, breaks=c(0, .001, .01, .05, .1), labels=c("***","**","*","^")),
             seurat_label_f= factor(seurat_label_f, levels=levels(tmp$x_labels))
      )
    tmp.stats = filter(plot_both, key_genes==plot_gene, annotation=="Seurat labels", condition %in% unlist(strsplit(target_group, "\\.")), sex==target_sex) %>%
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
    arrangeGrob(bp1, bp2, ncol=2, top=textGrob(plot_gene, gp=gpar(fontsize=14,fontface=4)))
  } else {
    arrangeGrob(bp1+ylim(ylim2), bp2+ylim(ylim2), ncol=2, top=textGrob(plot_gene, gp=gpar(fontsize=14,fontface=4)))
  } 
}
