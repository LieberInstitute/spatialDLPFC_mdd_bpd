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
deg.modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules.csv")

#pick out top predictor DEG network subset
mod_subset = c("SNHG14","PRKAR1A","COX4I1","EEF1A1",
  "PLP1",
  "GFAP","GLUL","HSPA1A","MT1M",
        "JUNB","ADAMTS1",
  "IFITM3","CD74","A2M",
  "GRIN1","CAMK2N1","GAD1")

deg.subset = filter(deg.modules, TF %in% mod_subset)

cat("\nPlot unrefined modules in top predictor DEG modules to highlight...\n")
plist = lapply(mod_subset, function(x) {
  df1 = filter(deg.subset, TF==x) %>% arrange(desc(importance)) %>%
    select(TF, target, importance) %>% mutate(is_original=T)
  df2 = filter(deg.subset, target %in% df1$target) %>% select(TF, target, importance)

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

## define custom functions
refineTargets <- function(seed_df) {
  ### isolate out the least important module
  not.min = filter(seed_df, is_min_imp==F)
  ### the rest
  is.min = filter(seed_df, is_min_imp==T)
  
  ### if for any of the modules the importance is >1.5x the least important module, keep
  tmp3 = left_join(select(not.min, TF, target, importance, rank_imp),
                   select(is.min, TF, target, importance),
                   by="target", suffix=c("_not.min","_is.min")) %>%
    mutate(rel_imp=importance_not.min/importance_is.min)
  
  notpass = filter(tmp3, rel_imp<1.5)
  ## if the most important relationship didn't pass, keep all module-target relationships for that target
  keep.all = filter(seed_df, target %in% filter(notpass, rank_imp==1)$target) %>% 
    select(TF, target, importance)
  
  
  yespass = filter(tmp3, rel_imp>1.5)
  #### sanity check
  #nrow(filter(yespass, target %in% keep.all$target)) #0, good
  
  ## if only 1 relationship passed, keep and remove from future
  keep.1 = group_by(yespass, target) %>% add_tally() %>% filter(n==1)
  #### sanity check that all relationships in this bucket are the most importance module-target pair
  #nrow(filter(keep.1, rank_imp!=1)) #0, good
  
  keep.complete = bind_rows(keep.all, select(keep.1, TF=TF_not.min, target, importance=importance_not.min))
  remaining.multi = filter(yespass, !target %in% keep.1$target) %>%
    group_by(target) %>% add_tally() %>% filter(n>1)
  
  return(list("keep"=keep.complete, "evaluate"=filter(not.min, target %in% remaining.multi$target)))
}


## isolate genes that are unique to 1 module
unique.targets = group_by(deg.subset, target) %>% add_tally(name="n_mods") %>% filter(n_mods==1) %>%
  select(TF, target, importance)

## isolate genes that are present in >1 module
multi.targets = group_by(deg.subset, target) %>% add_tally(name="n_mods") %>% filter(n_mods>1) %>%
  group_by(target) %>% mutate(rank_imp=rank(-importance), is_min_imp=rank_imp==n_mods)

## set unique as seed for keeping refined module-target pairs
refined.targets <- unique.targets
## while there are multi targets that have unequal importance between different module-target pairs, run refineTargets function
while(nrow(multi.targets)>0) {
  cat("\nDim refined.targets:", nrow(refined.targets))
  outList <- refineTargets(seed_df= multi.targets)
  refined.targets <- bind_rows(refined.targets, outList$keep)
  multi.targets <- group_by(outList$evaluate, target) %>% mutate(n_mods=n(), rank_imp=rank(-importance), is_min_imp=rank_imp==n_mods) 
}

## filter out module genes from targets (unfair because strongest importance should be its own module but Inf weight not evaluatable)
refined.targets <- filter(refined.targets, !target %in% mod_subset)

# plot refined modules
cat("\nPlot refined modules in top predictor DEG modules to highlight...\n")
plist2 = lapply(mod_subset, function(x) {
  df1 = filter(refined.targets, TF==x) %>% arrange(desc(importance)) %>%
    mutate(is_original=T)
  df2 = filter(refined.targets, target %in% df1$target)

  df2 = mutate(df2, x_lab= factor(target, levels=df1$target)) %>%
    left_join(df1, by=c("TF","target","importance")) %>%
    mutate(is_original= ifelse(is.na(is_original), "F", "T"))
  ggplot(df2, aes(x=x_lab, y=importance, color=is_original))+
    geom_point(size=.5)+scale_color_manual(values=c("F"="red","T"="black"))+
    scale_y_continuous(limits=c(0,max(df2$importance)))+labs(title=x, y="importance", x="module target genes")+
    theme_minimal()+theme(axis.text.x=element_blank(), legend.position="none")
})

## number of DEGs represented
source("code/09_DEG_GRN/load_DEGs.r")
mbv.degs = unique(sig.df$gene_name)
cat("\nNumber of DEGs present in the", length(mod_subset), "refined, un-filtered modules:", length(intersect(mbv.degs, c(refined.targets$target, refined.targets$TF))),"\n")

## evaluate module size, minimum size =10 
(mod.size = group_by(refined.targets, TF) %>% tally(name="n_targets") %>%
  arrange(n_targets))
too.small = filter(mod.size, n_targets<10)

if(nrow(too.small)>0) {
	cat("\n\nRemove modules with <10 DEG targets:", unique(too.small$TF),"\n")
	(mod_subset = setdiff(mod_subset, unique(too.small$TF)))
	
	cat("\n\nRe-run refining to distribute these targets to other modules (if link present)...\n")
	deg.subset = filter(deg.modules, TF %in% mod_subset)
	
	## isolate genes that are unique to 1 module
	unique.targets = group_by(deg.subset, target) %>% add_tally(name="n_mods") %>% filter(n_mods==1) %>%
	  select(TF, target, importance)

	## isolate genes that are present in >1 module
	multi.targets = group_by(deg.subset, target) %>% add_tally(name="n_mods") %>% filter(n_mods>1) %>%
	  group_by(target) %>% mutate(rank_imp=rank(-importance), is_min_imp=rank_imp==n_mods)

	## set unique as seed for keeping refined module-target pairs
	refined.targets <- unique.targets
	## while there are multi targets that have unequal importance between different module-target pairs, run refineTargets function
	while(nrow(multi.targets)>0) {
	  cat("\nDim refined.targets:", nrow(refined.targets))
	  outList <- refineTargets(seed_df= multi.targets)
	  refined.targets <- bind_rows(refined.targets, outList$keep)
	  multi.targets <- group_by(outList$evaluate, target) %>% mutate(n_mods=n(), rank_imp=rank(-importance), is_min_imp=rank_imp==n_mods)
	}

	## filter out module genes from targets (unfair because strongest importance should be its own module but Inf weight not evaluatable)
	refined.targets <- filter(refined.targets, !target %in% mod_subset)

	## number of DEGs represented
	cat("\nNumber of DEGs present in the", length(mod_subset), "refined, filtered modules:", length(intersect(mbv.degs, c(refined.targets$target, refined.targets$TF))),"\n")

	## look at new module sizes
	(mod.size = group_by(refined.targets, TF) %>% tally(name="n_targets") %>%
	  arrange(n_targets))
}

refined.modules <- left_join(refined.targets, deg.subset[,1:5], by=c("TF","target","importance"))
write.csv(refined.modules,
          "processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv", row.names=F)
cat("\nSaved refined top predictor DEG modules adj. file to: processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv\n")

# plot refined modules
cat("\nPlot refined modules in top predictor DEG modules to highlight...\n")
plist3 = lapply(mod_subset, function(x) {
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

saveRDS(list("igraph"=g, "layout"=lay1), "processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_DEG-modules-refined_igraph.rda")
cat("\nSaved refined DEG module igraph to: processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_DEG-modules-refined_igraph.rda\n")

pdf(file="plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_refining-DEG-modules.pdf",
    width=7, height=9)
do.call(grid.arrange, c(plist, ncol=3, top="Before refining..."))
do.call(grid.arrange, c(plist2, ncol=3, top="After refining..."))
do.call(grid.arrange, c(plist3, ncol=3, top="After refining and filtering..."))
plot(g, layout=lay1, main="Refined DEG modules")
dev.off()

cat("\nSaved all plots to: plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_refining-DEG-modules.pdf\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
