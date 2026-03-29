setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(igraph)
	library(gridExtra)
	library(fgsea)
})

set.seed(123)

# load in module sets
modules2 = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules.csv")

#pick out top predictor DEG network subset
mod_subset = c("UQCRH","PRKAR1A","CAMK2N1","GRIN1","GAD1",
               "SNHG14","GLUL","A2M","IFITM3","ADAMTS1","CD74",
               "FTL","PLP1","APLP1","HSPA1A","MT1X")

modules3 = filter(modules2, TF %in% mod_subset)

cat("\nPlot unrefined modules in top predictor DEG modules to highlight...\n")
plist = lapply(mod_subset, function(x) {
  df1 = filter(modules3, TF==x) %>% arrange(desc(importance)) %>%
    select(TF, target, importance) %>% mutate(is_original=T)
  df2 = filter(modules3, target %in% df1$target) %>% select(TF, target, importance)
  
  df2 = mutate(df2, x_lab= factor(target, levels=df1$target)) %>%
    left_join(df1, by=c("TF","target","importance")) %>% 
    mutate(is_original= ifelse(is.na(is_original), "F", "T")) 
  ggplot(df2, aes(x=x_lab, y=importance, color=is_original))+
    geom_point(size=.5)+scale_color_manual(values=c("F"="red","T"="black"))+
    scale_y_continuous(limits=c(0,max(df2$importance)))+labs(title=x, y="importance", x="module target genes")+
    theme_minimal()+theme(axis.text.x=element_blank(), legend.position="none")
})


#refine subset of DEG modules
cat("\nRefine top predictor DEG modules to reduce the amount of overlap of non-important target genes...\n")

unique.targets = group_by(modules3, target) %>% add_tally(name="n_mods") %>% filter(n_mods==1) %>%
  select(TF, target, importance, n_mods)

multi1 = group_by(modules3, target) %>% add_tally(name="n_mods") %>% filter(n_mods>1)
tmp = group_by(multi1, target) %>% summarise(max_imp=max(importance), min_imp=min(importance)) %>%
  arrange(desc(max_imp))

multi1 = left_join(multi1, tmp) %>% mutate(rel_min=importance/min_imp) %>%
  select(TF, target, importance, n_mods, rel_min) %>%
  mutate(x_lab=factor(target, levels=tmp$target))

## plot used to pick 1.5x min importance as cutoff
#ggplot(multi1, aes(x=x_lab, y=importance, color=cut(rel_min, breaks=c(0,1.5,2,100))))+
#  geom_point(size=.5)+
#  scale_color_manual(values=c("grey","red","black"))+
#  theme_minimal()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=8))

# if target has module where importance is >1.5x the min importance, only keep 
# target as module component in those cases
# if target doesn't have module where importance is >1.5 the min importance, keep all
filter.mins = group_by(multi1, target) %>% summarise(n_keep=sum(rel_min>1.5))
keep.all.mods = filter.mins$target[filter.mins$n_keep==0]
keep.top.mods = filter.mins$target[filter.mins$n_keep>0]

multi1 = bind_rows(filter(multi1, target %in% keep.all.mods),
                   filter(multi1, target %in% keep.top.mods, rel_min>1.5)) %>%
  select(TF, target, importance) %>%
  group_by(target) %>% add_tally(name="n_mods")

unique.targets = bind_rows(unique.targets, filter(multi1, n_mods==1))

multi2 = filter(multi1, n_mods>1)
tmp = group_by(multi2, target) %>% summarise(max_imp=max(importance), min_imp=min(importance)) %>%
  arrange(desc(max_imp))

multi2 = left_join(multi2, tmp) %>% mutate(rel_min=importance/min_imp) %>%
  select(TF, target, importance, n_mods, rel_min)

unique.targets = bind_rows(unique.targets, 
                           filter(multi2, rel_min>1.5) %>% mutate(n_mods=1) %>% 
                             select(TF, target, importance, n_mods))

multi.targets = filter(multi2, rel_min<= 1.5, !target %in% unique.targets$target) %>% 
  select(TF, target, importance, n_mods)


refined.modules = bind_rows(unique.targets, multi.targets) %>% filter(!target %in% mod_subset)

#group_by(refined.modules, target) %>% add_tally(name="check_n") %>%
#  filter(n_mods!=check_n)
## good, everything matches

# save refined modules
write.csv(left_join(refined.modules[,1:3], modules2[,1:5], by=c("TF","target","importance")),
          "processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv", row.names=F)
cat("\nSaved refined top predictor DEG modules adj. file to: processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv\n")


# plot refined modules
cat("\nPlot refined modules in top predictor DEG modules to highlight...\n")
plist2 = lapply(mod_subset, function(x) {
  df1 = filter(refined.modules, TF==x) %>% arrange(desc(importance)) %>%
    mutate(is_original=T)
  df2 = filter(refined.modules, target %in% df1$target) 

  df2 = mutate(df2, x_lab= factor(target, levels=df1$target)) %>%
    left_join(df1, by=c("TF","target","importance","n_mods")) %>% 
    mutate(is_original= ifelse(is.na(is_original), "F", "T")) 
  ggplot(df2, aes(x=x_lab, y=importance, color=is_original))+
    geom_point(size=.5)+scale_color_manual(values=c("F"="red","T"="black"))+
    scale_y_continuous(limits=c(0,max(df2$importance)))+labs(title=x, y="importance", x="module target genes")+
    theme_minimal()+theme(axis.text.x=element_blank(), legend.position="none")
})

# make refined igraph
cat("\nPlot refined modules igraph...\n")
e.df = select(ungroup(refined.modules), from=TF, to=target, weight=importance)
n.df = data.frame("node"=union(e.df$from, e.df$to),
	"is_TF"= union(e.df$from, e.df$to) %in% e.df$from)

g <- graph_from_data_frame(e.df, directed=T, vertices=n.df)
E(g)$weight = e.df$weight
E(g)$arrow.size = 0
E(g)$color = "gray"
  
V(g)$size = ifelse(n.df$is_TF, 5, 1)
V(g)$label = ifelse(n.df$is_TF, n.df$node, "")
V(g)$label.color = "black"
V(g)$color = "lightgoldenrod"
  
V(g)$label.family="sans"
V(g)$frame.width=0
V(g)$frame.color=NA
  
set.seed(123)
lay1= layout_with_kk(g, weights=sqrt(E(g)$weight))
set.seed(123)

saveRDS(list("igraph"=g, "layout"=lay1), "processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_top-predictor-DEG-modules-refined.rda")
cat("\nSaved refined DEG module igraph to: processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_top-predictor-DEG-modules-refined.rda\n")

# perform GSEA based on logFC rankings for enrichment of DEG modules
cat("\n\nGSEA on refined top predictor DEG modules to highlight...\n")
modList = mod_subset
names(modList) = mod_subset
modList = lapply(modList, function(x) c(x, filter(refined.modules, TF==x)$target))
sapply(modList, length)

runFGSEA <- function(filtered_df, .gmtList, adjp_threshold=.05, seed=123, ...) {
  set.seed(seed)
  #cat("\nSeed check:", seed.check(),"\n")
  lf.df = group_by(filtered_df, gene_name) %>% summarise(avg_lf=mean(logFC)) %>% 
    arrange(desc(avg_lf))
  lf.order = lf.df$avg_lf
  names(lf.order) = lf.df$gene_name
  results = fgseaMultilevel(.gmtList, stats=lf.order, scoreType="std", minSize=10, maxSize=500, ...)
  set.seed(seed)
  sig = filter(results, padj<adjp_threshold)
  sig$leadingEdge2 = sapply(sig$leadingEdge, paste, collapse="/")
  #cat("Seed check:", seed.check(),"\n")
  return(sig)
}


la.sm = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", 
                        "smoothed-k9-1663", "_rev-gene-input_moderated-t-test.csv"))
la.se = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", 
                        "seurat-pc30", "_rev-gene-input_moderated-t-test.csv"))

both.genes = intersect(la.sm$gene_id, la.se$gene_id)

la.both = full_join(filter(la.sm, gene_id %in% both.genes) %>% select(gene_id, gene_name, coef, logFC),
                    filter(la.se, gene_id %in% both.genes) %>% select(gene_id, gene_name, coef, logFC),
                    by=c("gene_name","gene_id","coef"), suffix=c("_sm","_se"))
la.both$logFC = (la.both$logFC_sm+la.both$logFC_se)/2

comparisons = c("F_NTC.MDD","F_NTC.BPD","F_MDD.BPD","M_NTC.BPD","M_MDD.BPD")
names(comparisons) <- comparisons

res1 = lapply(comparisons, function(x) {
  df1 = filter(la.both, coef==x)
  out1 = runFGSEA(df1, modList, adjp_threshold = 1)
  return(mutate(out1, sex.group=x))
})

#save results
saveRDS(res1, "processed-data/09_DEG_GRN/fgsea_refined-modules_LA-logFC-averaged.rda")
cat("\nSaved GSEA output to: processed-data/09_DEG_GRN/fgsea_refined-modules_LA-logFC-averaged.rda\n")

# print sig results
do.call(rbind, res1) %>% filter(padj<.05) %>%
  select(sex.group, pathway, NES) %>%
  tidyr::pivot_wider(names_from="sex.group", values_from="NES", values_fill=0) %>%
  mutate(pathway=factor(pathway, levels=mod_subset)) %>% arrange(pathway)

# plot enrichment as dotplot
plot.df = do.call(rbind, res1) %>% 
  mutate(is_sig= factor(padj<.05, levels=c("FALSE","TRUE"),
                        labels=c("NS","padj<.05")),
         y_lab=factor(pathway, levels=rev(mod_subset)),
         x_lab=factor(sex.group, levels=comparisons))
#max(abs(tmp$NES))

p1 <- ggplot(plot.df, aes(x=x_lab, y=y_lab, fill=NES, size=is_sig))+
  geom_count(shape=21, color="black")+scale_size_manual(values=c(2,6))+
  scale_fill_gradientn(colors=colorRampPalette(RColorBrewer::brewer.pal(n=7,"RdYlBu")[7:1])(100),
                        limits=c(-3.2,3.2))+
  geom_vline(aes(xintercept=3.5), lty=2)+
  labs(y="Top predictor DEG module", x="L-A logFC comparison")+
  theme(axis.text.y=element_text(face="italic", size=12))


cat("\nPlot GSEA enrichment curves...\n")
plotGSEA_jt <- function(dxsex, DEG_module) {
  r1 = filter(la.both, coef==dxsex)
  lf.df = group_by(r1, gene_name) %>% summarise(avg_lf=mean(logFC)) %>% 
    arrange(desc(avg_lf))
  lf.order = lf.df$avg_lf
  names(lf.order) = lf.df$gene_name
  ptitle = paste0(DEG_module, " (", dxsex, ")")
  is_sig = res1[[dxsex]]$padj[res1[[dxsex]]$pathway==DEG_module]<.05
  if(is_sig) {
    psub = factor(sign(res1[[dxsex]]$NES[res1[[dxsex]]$pathway==DEG_module]),
                  levels=c(-1,1), labels=c("Depleted","Enriched"))
  } else {
    psub = "Not significant"
  }
  p1 <- plotEnrichment(modList[[DEG_module]], lf.order)+labs(title=ptitle, subtitle=psub)
  return(p1)
}

lay_mat= cbind(c(1,2,3),c(NA,4,5))

pdf(file="plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_subset-refined-modules-and-GSEA.pdf", 
    width=7, height=9)
do.call(grid.arrange, c(plist, ncol=3, top="Before refining..."))
do.call(grid.arrange, c(plist2, ncol=3, top="After refining..."))
plot(g, layout=lay1, main="Refined DEG modules")
p1
for(i in setdiff(mod_subset,"ADAMTS1")) {
  plist3 = lapply(comparisons, plotGSEA_jt, DEG_module=i)
  plot(arrangeGrob(grobs=plist3, layout_matrix=lay_mat))
}
dev.off()
cat("\nSaved all plots to: plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_subset-refined-modules-and-GSEA.pdf\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
