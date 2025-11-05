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

plotViolin <- function(plot_gene, plot_both, DE_list, color_by=c("dx","cluster")) {
  adj.results_sm = DE_list[["LA_smoothed"]]
  adj.results_se = DE_list[["LA_seurat"]]
  
  if(color_by=="dx") plot_both$fill_color= plot_both$condition
  if(color_by=="cluster") plot_both$fill_color= plot_both$x_labels
  
  suppressMessages({
    tmp = filter(plot_both, key_genes==plot_gene)
    ylim1 = c(min(plot_both$logcounts), max(tmp$logcounts))
    
    #smoothed
    if(color_by=="dx") fillpal = cpList$dx.pal
    if(color_by=="cluster") fillpal = cpList$smoothed.bright
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
                 aes(x=condition, y=logcounts))+
      ggbeeswarm::geom_quasirandom(aes(color=fill_color))+
      scale_color_manual(values=fillpal)+
      geom_boxplot(fill="white", alpha=.5, outliers = F)+
      labs(title="PRECAST (smoothed)", subtitle=subt)+ylim(ylim1)+
      theme_bw()+theme(legend.position="none", legend.title = element_blank(), axis.title.x=element_blank(),
                       plot.subtitle=element_text(color=subt_color))
    
    #seurat labels
    if(color_by=="dx") fillpal = cpList$dx.pal
    if(color_by=="cluster") {
      fillpal = cpList$transfer.bright
      names(fillpal)[1] = "M.V"
    }
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
                 aes(x=condition, y=logcounts))+
      ggbeeswarm::geom_quasirandom(aes(color=fill_color))+
      scale_color_manual(values=fillpal)+
      geom_boxplot(fill="white", alpha=.5, outliers = F)+
      labs(title="Seurat labels", subtitle=subt)+ylim(ylim1)+
      theme_bw()+theme(legend.position="none", legend.title = element_blank(), axis.title.x=element_blank(),
                       plot.subtitle=element_text(color=subt_color))
    
  })
  
  arrangeGrob(bp1, bp2, ncol=2, top=textGrob(plot_gene, gp=gpar(fontsize=14,fontface=4)))
}
