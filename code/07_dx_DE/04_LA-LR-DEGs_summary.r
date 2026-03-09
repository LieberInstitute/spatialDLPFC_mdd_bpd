setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(scater)
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
	library(ggVennDiagram)
})
set.seed(123)

comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")
names(comparisons) <- comparisons

#results_set = "smoothed-k9-1663"
#load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
#spe_pseudo$cluster = spe_pseudo$smoothed_k9_1663
#comp_names = c("L1","L2","L3dot4","L5","L6","WM")
#names(comp_names) = c("L1","L2","L3.4","L5","L6","WM")

#results_set = "seurat-pc30"
#load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
#spe_pseudo$cluster = spe_pseudo$seurat_label
#comp_names = c("MicrodotVasc","Astro","L2dot3","L4","Inhb","L5","L6","Oligo")
#names(comp_names) = c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo")

results_set = "custom-cluster"
load("processed-data/06_pseudobulk/custom_cluster/spe_n119_pseudo_sample-custom-cluster_norm-filt.Rdata")
spe_pseudo$cluster = spe_pseudo$custom_cluster
comp_names = c("MicrodotVasc","AstrodotL1","AstrodotNrn","L2","L3","L4","Inhb","L5","L6","WM")
names(comp_names) = c("Micro.Vasc","Astro.L1","Astro.Nrn","L2","L3","L4","Inhb","L5","L6","WM")

# L-A ----
## f test
f.df = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", 
                       results_set, "_rev-gene-input_F-test.csv"))

f_sig = filter(f.df, adj.P.Val<.05)

# prop detected per capture area
load("processed-data/06_pseudobulk/spe_n119_pseudo-dotplot_sample-id.Rdata")

f_sig$med_prop.spots.detected = round(rowMedians(assay(spe_summ, "logcounts.prop.detected")[f_sig$gene_id,]),3)


## t test
t.df = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", 
                       results_set, "_rev-gene-input_moderated-t-test.csv"))

t_sig_pvt = filter(t.df, gene_id %in% f_sig$gene_id) %>% 
  group_by(gene_id, gene_name) %>% mutate(n_ttest_sig=sum(adj.P.Val<.05)) %>%
  ungroup() %>%
  #mutate(t_sig= factor(adj.P.Val<.05, levels=c("FALSE","TRUE"), labels=c("NS","padj<.05"))) %>% 
  mutate(t_sig= cut(adj.P.Val, breaks=c(0,.0001, .01, .05, 1), labels=c("padj<.0001","padj<.01","padj<.05","NS"))) %>%
  select(AveExpr,n_ttest_sig, t_sig, gene_id, gene_name, coef) %>%
  tidyr::pivot_wider(names_from="coef", values_from="t_sig")

la.degs = left_join(f_sig, t_sig_pvt, by=c("AveExpr","gene_id","gene_name"), suffix=c("_coef","_ttest"))

write.csv(la.degs, file=paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,"_dx-sex_degs-F-test-t-test.csv"), row.names=F)
cat("\nSaved L-A DEGs to:", paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,"_dx-sex_degs-F-test-t-test.csv"), "\n\n")

## plot
tmp = filter(t.df, gene_id %in% f_sig$gene_id, adj.P.Val<.05) %>% 
  mutate(coef=factor(coef, levels=comparisons),
         dir=factor(sign(logFC), levels=c(-1,1), labels=c("dn","up"))) %>%
  group_by(coef, dir, .drop=F) %>% tally()
p1 <- ggplot(tmp, aes(x=coef, y=n))+ 
  geom_bar(data=filter(tmp, dir=="up"), stat="identity", fill="#FFC0B5", color="tomato", width=.8)+
  geom_bar(data=filter(tmp, dir=="dn"), aes(y=-n), stat="identity", fill="#CFEBF7", color="skyblue", width=.8)+
  scale_x_discrete("", labels=gsub("\\.", " vs. ", gsub("_","\n",comparisons)))+
  labs(y="# genes (F test adj. p<.05 t-test adj. p<.05)", title=paste0("L-A model: ", results_set),
       subtitle=paste0("Of ", nrow(la.degs), " genes with F test adj. p<.05, ", 
                       sum(la.degs$n_ttest_sig>0), " genes also with t-test adj. p<.05."))+
  theme_minimal()


# L-R ----
f.df_lr = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", 
                       results_set, "_rev-gene-input_F-test_all.csv"))
#remove cluster="all" column that isn't helpful for just this dframe and could be confusing in the future
f.df_lr = f.df_lr[,setdiff(colnames(f.df_lr),"cluster")]

f_sig_lr = filter(f.df_lr, adj.P.Val<.05)

## prop detected
if(results_set=="smoothed-k9-1663") {
  load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo-dotplot_sample-id-smoothed-n1663-k9.Rdata")
  spe_summ$cluster = spe_summ$smoothed_k9_1663
}
if(results_set=="seurat-pc30") {
  load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-dotplot_sample-id-seurat-pc30.Rdata")
  spe_summ$cluster = spe_summ$seurat_label
}
if(results_set=="custom-cluster") {
  load("processed-data/06_pseudobulk/custom_cluster/spe_n119_pseudo-dotplot_sample-id-custom-cluster.Rdata")
  spe_summ$cluster = spe_summ$custom_cluster
}

df = do.call(rbind, lapply(names(comp_names), function(x) {
  rmeds = rowMedians(assay(spe_summ, "logcounts.prop.detected")[f_sig_lr$gene_id,spe_summ$cluster==x])
  data.frame("gene_id"=f_sig_lr$gene_id, "cluster"=x, med_prop.spots.detected=round(rmeds,3))
}))
prop.df = tidyr::pivot_wider(df, names_from="cluster", names_prefix="med_prop.spots.detected_", values_from="med_prop.spots.detected")

f_sig_lr = left_join(f_sig_lr, prop.df)

# ttest p values
t.df_lr = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", 
                       results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
  mutate(cluster=factor(cluster, levels=names(comp_names)),
         sex.group=factor(paste(sex, group, sep="_"), levels=comparisons))

t_sig_pvt_lr = filter(t.df_lr, gene_id %in% f_sig_lr$gene_id) %>% 
  group_by(gene_id, gene_name) %>% mutate(n_ttest_sig=sum(adj.P.Val<.05)) %>%
  ungroup() %>%
#  mutate(t_sig= factor(adj.P.Val<.05, levels=c("FALSE","TRUE"), labels=c("NS","padj<.05"))) %>% 
  mutate(t_sig= cut(adj.P.Val, breaks=c(0,.0001, .01, .05, 1), labels=c("padj<.0001","padj<.01","padj<.05","NS"))) %>%
  select(AveExpr,n_ttest_sig, t_sig, gene_id, gene_name, coef) %>%
  tidyr::pivot_wider(names_from="coef", values_from="t_sig")

sig.by.cluster = filter(t.df_lr, gene_id %in% f_sig_lr$gene_id) %>% 
  group_by(gene_id, gene_name, cluster) %>% summarise(n_ttest_sig=sum(adj.P.Val<.05)) %>%
  tidyr::pivot_wider(names_from="cluster", names_prefix="n_ttest_sig_", values_from="n_ttest_sig")
sig.by.group = filter(t.df_lr, gene_id %in% f_sig_lr$gene_id) %>% 
  group_by(gene_id, gene_name, sex.group) %>% summarise(n_ttest_sig=sum(adj.P.Val<.05)) %>%
  tidyr::pivot_wider(names_from="sex.group", names_prefix="n_ttest_sig_", values_from="n_ttest_sig")

t_sig_pvt_lr2 = left_join(t_sig_pvt_lr, sig.by.cluster) %>% left_join(sig.by.group)
t_sig_pvt_lr2 = t_sig_pvt_lr2[,c("AveExpr","gene_id","gene_name","n_ttest_sig",
                                 colnames(sig.by.cluster)[3:ncol(sig.by.cluster)],
                                 colnames(sig.by.group)[3:ncol(sig.by.group)],
                                 colnames(t_sig_pvt_lr)[5:ncol(t_sig_pvt_lr)])]

lr.degs = left_join(f_sig_lr, t_sig_pvt_lr2, by=c("AveExpr","gene_id","gene_name"), suffix=c("_coef","_ttest"))

write.csv(lr.degs, file=paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,"_dx-sex_degs-F-test-t-test.csv"), row.names=F)
cat("\nSaved L-R DEGs to:", paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,"_dx-sex_degs-F-test-t-test.csv"), "\n\n")


## venn diagram of overlap
vp <- ggVennDiagram(list("LA"=la.degs$gene_id, "LR"=lr.degs$gene_id), label="count",
                    edge_size=.5, set_color=c("black","darkgreen"),
              category.names= c("L-A genes (F adj. p<.05)", "L-R genes (F adj. p<.05)"))+
  scale_fill_gradient(low="white", high="grey60")+
  labs(title=results_set)+
  xlim(-5,5)+
  coord_flip()+theme(legend.position="none", plot.title=element_text(hjust=.5, color="grey"),
                     plot.margin=margin(6,50,10,50, "pt"), aspect.ratio=.8,
                     text=element_text(size=10))

b <- ggplot_build(vp)  
b$data[[3]][,"x"] = c(5,-5)
b$data[[3]][,"y"] = c(0,4)

## plot
tmp = filter(t.df_lr, gene_id %in% f_sig_lr$gene_id, adj.P.Val<.05) %>% 
  mutate(dir=factor(sign(logFC), levels=c(-1,1), labels=c("dn","up"))) %>%
  group_by(sex.group, cluster, dir, .drop=F) %>% tally()
p2 <- ggplot(tmp, aes(x=cluster, y=n))+ 
  geom_bar(data=filter(tmp, dir=="up"), stat="identity", fill="#FFC0B5", color="tomato")+
  geom_bar(data=filter(tmp, dir=="dn"), aes(y=-n), stat="identity", fill="#CFEBF7", color="skyblue")+
  facet_wrap(vars(sex.group))+
  labs(y="# genes (F test adj. p<.05 t-test adj. p<.05)", 
       title=paste0("L-R model: ", results_set), x="",
       subtitle=paste0("Of ", nrow(lr.degs), " genes with F test adj. p<.05, ", 
                       sum(lr.degs$n_ttest_sig>0), " genes also with t-test adj. p<.05."))+
  theme_minimal()+
  theme(panel.border = element_rect(fill=NA, color="grey80"),
        strip.background = element_rect(color="grey80", fill="grey80"))

if(results_set=="seurat-pc30") p2 <- p2+scale_x_discrete(labels=c("M.V","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))
if(results_set=="custom-cluster") p2 <- p2+scale_x_discrete(labels=c("M.V","Ast\nL1","Ast\nNrn","L2","L3","L4","Inhb","L5","L6","WM"))

## plot each layers unique genes
source("code/06_pseudobulk/custom_functions.r")
cond_sex = c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M")

if(results_set=="smoothed-k9-1663") {
  load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo-dotplot_dx-sex-smoothed-n1663-k9.Rdata")
  spe_summ$cluster = factor(spe_summ$smoothed_k9_1663, levels=names(comp_names))
}
if(results_set=="seurat-pc30") {
  load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-dotplot_dx-sex-seurat-pc30.Rdata")
  spe_summ$cluster = factor(spe_summ$seurat_label, levels=names(comp_names))
}
if(results_set=="custom-cluster") {
  load("processed-data/06_pseudobulk/custom_cluster/spe_n119_pseudo-dotplot_dx-sex-custom-cluster.Rdata")
  spe_summ$cluster = factor(spe_summ$custom_cluster, levels=names(comp_names))
}

spe_summ$sample_id = factor(paste(spe_summ$condition, spe_summ$sex, spe_summ$cluster),
                            levels=as.character(outer(cond_sex, names(comp_names), paste)))
colnames(spe_summ) <- spe_summ$sample_id

#set color limits based on all of the LR degs being plotted to color scale is consistent across plots
lim.df = dotplotDF(spe_summ, lr.degs$gene_name[lr.degs$n_ttest_sig>0], summarize_groups=T, swap_rownames="gene_name",
                 cluster_labels="cluster", row_data=NULL)
color_limits = round(max(abs(lim.df$mean_expr_scaled)),1)
if(color_limits<max(abs(lim.df$mean_expr_scaled))) color_limits = color_limits+.1

plist = lapply(names(comp_names), function(x) {
  # select layer DEGs
  plot.genes = lr.degs$gene_name[lr.degs[[paste0("n_ttest_sig_", x)]]>0]
  
  # determine order with heatmap
  hmp = plotGroupedHeatmap(spe_pseudo, features=plot.genes, swap_rownames="gene_name",
                           group="cluster", cluster_cols=F, silent=T, center=T)
  gene_order = rev(hmp$tree_row$label[hmp$tree_row$order])
  
  # if the number of layer DEGs is too long, split it into smaller sets for multiple plots
  if(length(plot.genes)>45) {
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
    
    # if splitting into two if still too large, split further
    if(other.clus>45) {
      hmp2 = plotGroupedHeatmap(spe_pseudo, features=other.set, swap_rownames="gene_name",
                               group="cluster", cluster_cols=F, silent=T, center=T)
      gene_order2 = rev(hmp2$tree_row$label[hmp2$tree_row$order])
      largest.group = other.clus
      .k=1
      while(largest.group>45) {
        .k=.k+1
        gene.grps = cutree(hmp2$tree_row, k=.k)
        largest.group = max(table(gene.grps))
      }
      mid.clus = table(gene.grps)[table(gene.grps)==largest.group]
      mid.set = names(gene.grps)[gene.grps==names(mid.clus)]
      splitGenes = c(splitGenes, list(intersect(gene_order2, mid.set)))
      other.clus = sum(table(gene.grps)[table(gene.grps)!=largest.group])
      if(other.clus>45) stop("Need to go further with splitting of dotplots.")
      other.set = intersect(gene_order2, names(gene.grps)[gene.grps!=names(mid.clus)])
    }
    splitGenes = c(splitGenes, list(other.set))
  } else {
    splitGenes = list(gene_order)
  }
  
  # now that the genes that are part of the dotplot are selected, make dotplots
  dotlist <- lapply(splitGenes, function(y) {
    #dotplot frame
    p.df = dotplotDF(spe_summ, y, summarize_groups=T, swap_rownames="gene_name",
                     cluster_labels="cluster", row_data=NULL) 
    p.df$gene_name_f = factor(p.df$gene_name, levels=y)
    #color_limits = round(max(abs(p.df$mean_expr_scaled)),1)
    #if(color_limits<max(abs(p.df$mean_expr_scaled))) color_limits = color_limits+.1

    #additional DF to label LR sig
    sig.join = mutate(t.df_lr, dir=factor(sign(logFC), levels=c(1,-1), labels=c("Inc.","Dec."))) %>%
	filter(cluster==x, adj.P.Val<.05, gene_name %in% y, sex.group %in% c("F_NTC.MDD","M_NTC.MDD","F_NTC.BPD","M_NTC.BPD")) %>%
	select(gene_name, gene_id, dir, cluster)
    check.n = distinct(sig.join, gene_id, gene_name, dir) %>% group_by(gene_id, gene_name) %>% tally()
    if(max(check.n$n)!=1) {
      both.dir = filter(check.n, n>1)
      sig.join.supp = sig.join[sig.join$gene_id %in% both.dir$gene_id,]
      sig.join.supp = sig.join.supp[1,]
      sig.join.supp[,"dir"] = "Both"
      sig.join <- bind_rows(filter(sig.join, !gene_id %in% both.dir$gene_id),
                            sig.join.supp)
    }
    sig.plot = right_join(p.df, sig.join, by=c("clusters"="cluster", "gene_name", "gene_id"))

    #generate plot
    p <- ggplot(p.df, aes(x=clusters, y=gene_name_f))+ 
                          #color=mean_expr_scaled, size=prop_spots))+
      geom_count(aes(fill=mean_expr_scaled, size=prop_spots), shape=21, color="grey")+
      scale_fill_gradient2(low="white", mid="lightgoldenrod", high="black",
                            midpoint=0, limits=c(-color_limits,color_limits))+
      scale_size(range=c(1,6), limits=c(0,1), breaks=c(0,.5,1))+
      geom_count(data=filter(sig.plot, dir=="Inc."), aes(size=0), color="tomato", show.legend = F)+
      geom_count(data=filter(sig.plot, dir=="Dec."), aes(size=0), color="dodgerblue", show.legend = F)+
      labs(fill="Avg. expr.\n(scaled)", size="Prop. of\nspots",
           y=paste0("F test adj. p<.05 & any ", x, " t-test adj. p<.05"), 
           title=paste0("L-R ", results_set, ": ", x, " DEGs"),
	   subtitle="Central red/blue dot indicates dir of DE vs NTC")+
      theme_minimal()+
      theme(axis.title.x=element_blank(), 
            #panel.background = element_rect(fill="grey90"),
            #panel.grid = element_line(color="grey80"),
            axis.text.y=element_text(face="italic"), legend.key.size=unit(15,"pt"),
            axis.title.y=element_text(margin=margin(0,20,0,20,"pt")))
    if(max(check.n$n)!=1) {
      p <- p+geom_count(data=filter(sig.plot, dir=="Both"), aes(size=0), color="purple", show.legend = F)+
        labs(subtitle="Central red/blue dot indicates dir of DE vs NTC (purple means sig each dir)")
    }
    return(p)
  })
  return(dotlist)
})


pdf(file=paste0("plots/07_dx_DE/dx-sex_pc3-age-nspots_", results_set,"_DEG-plots.pdf"))
grid.arrange(p1+theme(plot.margin=margin(6,50,6,50, "pt")), ggplot_gtable(b), ncol=1)
p2
for(i in 1:length(plist)) {
  if(length(plist[[i]])>1) {
    print(marrangeGrob(grobs=plist[[i]], ncol=1, nrow=1))
  } else {
    grid.arrange(plist[[i]][[1]])
  }
}
dev.off()
cat("\n\nPlots saved to:", paste0("plots/07_dx_DE/dx-sex_pc3-age-nspots_", results_set,"_DEG-plots.pdf"),"\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
