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
mod_subset = c("EEF1A1",#"EIF1",
	"PLP1","CD74","COX4I1","GLUL","IFITM3","GRIN1",
        "SNHG14","CAMK2N1","GAD1","A2M","PRKAR1A", #"UQCRH",
	"ADAMTS1","HSPA1A","MT1M","JUNB",
        "GFAP")

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

## isolate genes that are unique to 1 module
unique.targets = group_by(modules3, target) %>% add_tally(name="n_mods") %>% filter(n_mods==1) %>%
  select(TF, target, importance)

## isolate genes that are present in >1 module
multi.targets = group_by(modules3, target) %>% add_tally(name="n_mods") %>% filter(n_mods>1) %>%
	group_by(target) %>% mutate(rank_imp=rank(-importance), is_min_imp=rank_imp==n_mods)

## create first filter by isolating genes where the most important module is >1.5x more important than the second most important module
top.targets = filter(multi.targets, rank_imp==1)
second.targets = filter(multi.targets, rank_imp==2)

tmp = left_join(select(top.targets, TF, target, importance),
                select(second.targets, TF, target, importance, is_min_imp),
                 by=c("target"), suffix=c("_max","_min")) %>%
  mutate(rel_imp=importance_max/importance_min)

### if true, keep top only
keep.top = filter(tmp, rel_imp>1.5) %>% select(TF=TF_max, target, importance=importance_max)

### if not true AND there are only two modules for the gene, keep both
keep.both = bind_rows(filter(tmp, rel_imp<=1.5, is_min_imp==T) %>% select(TF=TF_max, target, importance=importance_max),
                     filter(tmp, rel_imp<=1.5, is_min_imp==T) %>% select(TF=TF_min, target, importance=importance_min))

## start dframe with refined module-targets
keep.complete = bind_rows(unique.targets, keep.top, keep.both)

## reduce the dframe with targets in >1 module to exclude elements we decided to keep 
multi.targets2 = filter(multi.targets, !target %in% union(keep.top$target, keep.both$target))

## create second filter by isolating genes where the most important module is >1.5x more important than the third most important module
top.targets = filter(multi.targets, rank_imp==1)
third.targets = filter(multi.targets2, rank_imp==3)

tmp1 = left_join(select(top.targets, TF, target, importance),
                select(third.targets, TF, target, importance, is_min_imp),
                by=c("target"), suffix=c("_max","_min")) %>%
  mutate(rel_imp=importance_max/importance_min)

### where the gene is present in only 3 modules and the most important is <1.5x more important than the least important, keep as target in all modules
keep.all = filter(multi.targets2, target %in% filter(tmp1, rel_imp<=1.5, is_min_imp==T)$target) %>% 
  select(TF, target, importance)

### where the most important module is >1.5x more important than the third most important module, keep the most important module
keep.top1 = filter(tmp1, rel_imp>1.5) %>% select(TF=TF_max, target, importance=importance_max)
#### in these cases, check to see if the second most important module is also >1.5x more important than the third most important module
tmp2 = left_join(select(third.targets, TF, target, importance, is_min_imp), 
          select(second.targets, TF, target, importance), by="target", suffix=c("_min","_max")) %>%
  mutate(rel_imp=importance_max/importance_min) %>%
  filter(target %in% keep.top1$target)
#### if second is also >1.5x, keep the second most important module
keep.top2 = filter(tmp2, rel_imp>1.5) %>% select(TF=TF_max, target, importance=importance_max)

## add kept decisions to all kept elements
keep.complete = bind_rows(keep.complete, keep.all, keep.top1, keep.top2)

#setdiff(multi.targets2$target, keep.complete$target)
#just a couple more to go that have rel_imp<1.5 and aren't yet min_imp

## there are a few targets remaining with 4-5 modules for each target, and we have already established that top 1-2 aren't greater than 3rd place
remaining.targets = filter(tmp1, rel_imp<1.5, is_min_imp==F)$target

### isolate out the least important module
not.min = filter(multi.targets2, target %in% remaining.targets, is_min_imp==F)
### the rest
is.min = filter(multi.targets2, target %in% remaining.targets, is_min_imp==T)

### if for any of the modules the importance is >1.5x the least important module, keep
tmp3 = left_join(select(not.min, TF, target, importance, rank_imp),
          select(is.min, TF, target, importance),
          by="target", suffix=c("_not.min","_is.min")) %>%
  mutate(rel_imp=importance_not.min/importance_is.min)

keep.any = filter(tmp3, rel_imp>1.5)

### if none are, keep all
keep.all = filter(multi.targets2, target %in% setdiff(remaining.targets, keep.any$target)) %>%
  select(TF, target, importance)

## final additions to refined set
keep.complete = bind_rows(keep.complete, select(keep.any, TF=TF_not.min, target, importance=importance_not.min),
          keep.all)

## remove any genes corresponding to our module subset from all other modules
### do this because weight is unfair, the weight/importance should always be highest for self but can't reflect that so just removing from other modules
refined.modules = filter(keep.complete, !target %in% mod_subset)

(mod.size = group_by(refined.modules, TF) %>% tally(name="n_targets") %>%
  arrange(n_targets))

too.small = filter(mod.size, n_targets<10)
refined.modules2 = filter(refined.modules, !TF %in% too.small$TF) %>%
  group_by(target) %>% add_tally(name="n_mods")

write.csv(left_join(refined.modules2[,1:3], modules2[,1:5], by=c("TF","target","importance")),
          "processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv", row.names=F)
cat("\nSaved refined top predictor DEG modules adj. file to: processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv\n")


# plot refined modules
cat("\nPlot refined modules in top predictor DEG modules to highlight...\n")
plist2 = lapply(mod_subset, function(x) {
  df1 = filter(refined.modules, TF==x) %>% arrange(desc(importance)) %>%
    mutate(is_original=T)
  df2 = filter(refined.modules, target %in% df1$target)

  df2 = mutate(df2, x_lab= factor(target, levels=df1$target)) %>%
    left_join(df1, by=c("TF","target","importance")) %>% 
    mutate(is_original= ifelse(is.na(is_original), "F", "T"))
  ggplot(df2, aes(x=x_lab, y=importance, color=is_original))+
    geom_point(size=.5)+scale_color_manual(values=c("F"="red","T"="black"))+
    scale_y_continuous(limits=c(0,max(df2$importance)))+labs(title=x, y="importance", x="module target genes")+
    theme_minimal()+theme(axis.text.x=element_blank(), legend.position="none")
})

# make refined igraph
cat("\nPlot refined modules igraph...\n")
e.df = select(ungroup(refined.modules2), from=TF, to=target, weight=importance)
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

saveRDS(list("igraph"=g, "layout"=lay1), "processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_DEG-modules-refined_igraph.rda")
cat("\nSaved refined DEG module igraph to: processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_DEG-modules-refined_igraph.rda\n")

pdf(file="plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_refining-DEG-modules.pdf",
    width=7, height=9)
do.call(grid.arrange, c(plist, ncol=3, top="Before refining..."))
do.call(grid.arrange, c(plist2, ncol=3, top="After refining..."))
plot(g, layout=lay1, main="Refined DEG modules")
dev.off()

cat("\nSaved all plots to: plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_refining-DEG-modules.pdf\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
