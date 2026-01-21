setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(SpatialExperiment)
	library(scater)
	library(gridExtra)
})
set.seed(123)

cpList <- readRDS("plots/colorPalettes.rds")
source("code/06_pseudobulk/custom_functions.r")

#results_set = "smoothed-k9-1663"
#comp_names = c("L1","L2","L3dot4","L5","L6","WM")
#names(comp_names) = c("L1","L2","L3.4","L5","L6","WM")
#col.pal = cpList$smoothed.bright

results_set = "seurat-pc30"
comp_names = c("MicrodotVasc","Astro","L2dot3","L4","Inhb","L5","L6","Oligo")
names(comp_names) = c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo")
col.pal = cpList$transfer.bright


comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")
names(comparisons) <- comparisons

cat("\n\n********************", results_set, "********************\n")

cat("\n> Layer-adjusted\n")
la.degs = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))
cat("\nNumber of L-A F test adj. p<.05:", nrow(la.degs), "\n")
cat("Number of L-A F test adj. p<.05 with sig t-test result:", sum(la.degs$n_ttest_sig>0), "\n")
cat("Number of L-A DEGs per dx*sex group:\n")
colSums(la.degs[,18:23]!="NS")

cat("\n\n> Layer-restricted\n")
lr.degs = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))
cat("\nNumber of L-R F test adj. p<.05:", nrow(lr.degs), "\n")
cat("Number of L-R F test adj. p<.05 with sig t-test result:", sum(lr.degs$n_ttest_sig>0), "\n")

tmp = colnames(lr.degs)[grep("n_ttest_sig", colnames(lr.degs))]
groups = tmp[grep("F_NTC.MDD", tmp):grep("M_MDD.BPD", tmp)]
cat("Number of L-R DEGs per dx*sex group:\n")
colSums(lr.degs[,groups]>0)
clusters = tmp[2:(grep("F_NTC.MDD", tmp)-1)]
cat("\nNumber of L-R DEGs per cluster:\n")
colSums(lr.degs[,clusters]>0)

cat("\n\n\n******* F test adj. p<.05 overlap only *******\n\n")
lr.and.la = intersect(la.degs$gene_name, lr.degs$gene_name)
cat("L-A and L-R overlap:", length(lr.and.la), "\n")

la.not.lr = setdiff(la.degs$gene_name, lr.degs$gene_name)
cat("L-A only:", length(la.not.lr), "\n")

lr.not.la = setdiff(lr.degs$gene_name, la.degs$gene_name)
cat("L-R only:", length(lr.not.la), "\n")
#scater::plotGroupedHeatmap(spe_save, features=lr.not.la_sm, group="smoothed_k9_1663", 
#                           center=T, show_rownames=F, cluster_cols=F, angle_col=0)

cat("\n\n\n******* DEGs (F test adj. p<.05 and t-test adj. p<.05) overlap *******\n\n")
la.degs_t = filter(la.degs, n_ttest_sig>0)
lr.degs_t = filter(lr.degs, n_ttest_sig>0)

lr.and.la_t = intersect(la.degs_t$gene_name, lr.degs_t$gene_name)
cat("L-A and L-R overlap:", length(lr.and.la_t), "\n")

la.not.lr_t = setdiff(la.degs_t$gene_name, lr.degs_t$gene_name)
cat("L-A only:", length(la.not.lr_t), "\n")

lr.not.la_t = setdiff(lr.degs_t$gene_name, la.degs_t$gene_name)
cat("L-R only:", length(lr.not.la_t),"\n")


cat("\nL-A and L-R overlap split by dx*sex group:\n")
sapply(comparisons, function(i) {
  any.clus = lr.degs_t[[paste0("n_ttest_sig_",i)]]>0
  la.res = la.degs_t[[paste0(i,"_ttest")]]!="NS"
  c("L-A_L-R_overlap"=length(intersect(lr.degs_t$gene_name[any.clus], la.degs_t$gene_name[la.res])),
    "LA_only"=length(setdiff(la.degs_t$gene_name[la.res], lr.degs_t$gene_name[any.clus])),
    "LR_only"=length(setdiff(lr.degs_t$gene_name[any.clus], la.degs_t$gene_name[la.res])))
})


# plot L-R overlap bar graph
lr.unique = lapply(comparisons, function(i) {
  any.clus = lr.degs_t[[paste0("n_ttest_sig_",i)]]>0
  la.res = la.degs_t[[paste0(i,"_ttest")]]!="NS"
  list("LR_only"=setdiff(lr.degs_t$gene_name[any.clus], la.degs_t$gene_name[la.res]),
       "LR_LA_overlap"=intersect(lr.degs_t$gene_name[any.clus], la.degs_t$gene_name[la.res]))
})

lrt = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_",
                      results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
  mutate(dir=factor(sign(logFC), levels=c(1,-1), labels=c("Inc.","Dec.")), 
         cluster=factor(cluster, levels=names(comp_names)),
         sex.group = factor(paste(sex, group, sep="_"), levels=comparisons)) 


plist <- lapply(comparisons, function(i) {
  tmp = filter(lrt, sex.group==i, adj.P.Val<.05)
  tmp.lr.overlap = filter(tmp, gene_name %in% lr.unique[[i]]$LR_LA_overlap) %>%
    group_by(cluster, .drop=F) %>% tally() %>% mutate(overlap="L-A & L-R overlap")
  tmp.lr.only = filter(tmp, gene_name %in% lr.unique[[i]]$LR_only) %>% 
    group_by(cluster) %>% tally() %>% mutate(overlap="L-R only")
  
  if(results_set=="smoothed-k9-1663") ymax=55
  if(results_set=="seurat-pc30") ymax=65
  
  p <- ggplot(bind_rows(tmp.lr.overlap, tmp.lr.only) %>% mutate(overlap=factor(overlap, levels=c("L-A & L-R overlap","L-R only"))),
         aes(x=cluster, y=n, fill=cluster))+
    geom_bar(stat="identity")+
    ylim(0,ymax)+facet_wrap(vars(overlap))+
    scale_fill_manual(values=col.pal, guide="none")+
    labs(y="# DEGs", x="", title=i)+
    theme_minimal()+theme(panel.border = element_rect(fill=NA, color="grey"),
                          plot.title=element_text(hjust=.5))
  if(results_set=="seurat-pc30") p <- p+scale_x_discrete(labels=c("M.V","A","L2.3","L4","In","L5","L6","Olg"))
  return(p)
})

# plot dotplot of uniquely L-R genes
cond_sex = c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M")

if(results_set=="smoothed-k9-1663") {
  load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
  spe_pseudo$cluster = spe_pseudo$smoothed_k9_1663
  load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo-dotplot_dx-sex-smoothed-n1663-k9.Rdata")
  spe_summ$cluster = factor(spe_summ$smoothed_k9_1663, levels=names(comp_names))
}
if(results_set=="seurat-pc30") {
  load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
  spe_pseudo$cluster = spe_pseudo$seurat_label
  load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-dotplot_dx-sex-seurat-pc30.Rdata")
  spe_summ$cluster = factor(spe_summ$seurat_label, levels=names(comp_names))
}

spe_summ$sample_id = factor(paste(spe_summ$condition, spe_summ$sex, spe_summ$cluster),
                            levels=as.character(outer(cond_sex, names(comp_names), paste)))
colnames(spe_summ) <- spe_summ$sample_id

plist2 <- lapply(comparisons, function(i) {
  plot.genes = lr.unique[[i]]$LR_only
  # determine order with heatmap
  hmp = plotGroupedHeatmap(spe_pseudo, features=plot.genes, swap_rownames="gene_name",
                           group="cluster", cluster_cols=F, silent=T, center=T,
                           main=i)
  gene_order = rev(hmp$tree_row$label[hmp$tree_row$order])
  # if the number of layer DEGs is too long, split it into smaller sets for multiple plots
  if(length(plot.genes)>50) {
    largest.group = length(plot.genes)
    .k=1
    while(largest.group>45) {
      .k=.k+1
      gene.grps = cutree(hmp$tree_row, k=.k)
      largest.group = max(table(gene.grps))
    }
    large.clus = table(gene.grps)[table(gene.grps)==largest.group]
    large.set = names(gene.grps)[gene.grps==names(large.clus)]
    
    splitGenes = list(intersect(gene_order, large.set))
    
    other.clus = sum(table(gene.grps)[table(gene.grps)!=largest.group])
    other.set = intersect(gene_order, names(gene.grps)[gene.grps!=names(large.clus)])
    splitGenes = c(splitGenes, list(other.set))
  } else {
    splitGenes = list(gene_order)
  }
  dotlist <- lapply(splitGenes, function(y) {
    #dotplot frame
    p.df = dotplotDF(spe_summ, y, summarize_groups=T, swap_rownames="gene_name",
                     cluster_labels="cluster", row_data=NULL) 
    p.df$gene_name_f = factor(p.df$gene_name, levels=y)
    color_limits = round(max(abs(p.df$mean_expr_scaled)),1)
    if(color_limits<max(abs(p.df$mean_expr_scaled))) color_limits = color_limits+.1
    
    #additional DF to label L-R sig
    sig.join = filter(lrt, sex.group==i, adj.P.Val<.05, gene_name %in% y) %>%
      select(gene_name, gene_id, dir, cluster)
    
    sig.plot = right_join(p.df, sig.join, by=c("clusters"="cluster", "gene_name", "gene_id"))
    
    #generate plot
    p <- ggplot(p.df, aes(x=clusters, y=gene_name_f))+
      geom_count(aes(color=mean_expr_scaled, size=prop_spots))+
      scale_color_gradient2(low="white", mid="lightgoldenrod", high="black",
                            midpoint=0, limits=c(-color_limits,color_limits))+
      scale_size(range=c(1,6), limits=c(0,1), breaks=c(0,.5,1))+
      geom_count(data=filter(sig.plot, dir=="Inc."), aes(size=0), color="tomato", show.legend = F)+
      geom_count(data=filter(sig.plot, dir=="Dec."), aes(size=0), color="dodgerblue", show.legend = F)+
      labs(color="Avg. expr.\n(scaled)", size="Prop. of\nspots",
           y="F test adj. p<.05 & t-test adj. p<.05", 
           title=paste0(results_set, " ", i, ": L-R only DEGs"))+
      theme_minimal()+
      theme(axis.title.x=element_blank(), 
            panel.background = element_rect(fill="grey90"),
            panel.grid = element_line(color="grey80"),
            axis.text.y=element_text(face="italic"), legend.key.size=unit(15,"pt"),
            axis.title.y=element_text(margin=margin(0,20,0,20,"pt")))
    return(p)
  })
  return(dotlist)
})

pdf(file=paste0("plots/07_dx_DE/dx-sex_pc3-age-nspots_", results_set, "_LA-LR-DEG-overlap.pdf"))
do.call(grid.arrange, c(plist, ncol=2))
for(i in 1:length(plist2)) {
  if(length(plist2[[i]])>1) {
    print(marrangeGrob(grobs=plist2[[i]], ncol=1, nrow=1))
  } else {
    grid.arrange(plist2[[i]][[1]])
  }
}
dev.off()

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
