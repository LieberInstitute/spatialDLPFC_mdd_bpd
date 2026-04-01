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

# remove modules with <10 targets
group_by(refined.modules, TF) %>% tally(name="n_targets") %>%
arrange(n_targets)

refined.modules2 = filter(refined.modules, TF!="ADAMTS1")

# save refined modules
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
    left_join(df1, by=c("TF","target","importance","n_mods")) %>% 
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

saveRDS(list("igraph"=g, "layout"=lay1), "processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_top-predictor-DEG-modules-refined.rda")
cat("\nSaved refined DEG module igraph to: processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_top-predictor-DEG-modules-refined.rda\n")


pdf(file="plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_subset-refined-modules.pdf",
    width=7, height=9)
do.call(grid.arrange, c(plist, ncol=3, top="Before refining..."))
do.call(grid.arrange, c(plist2, ncol=3, top="After refining..."))
plot(g, layout=lay1, main="Refined DEG modules")
dev.off()

cat("\nSaved all plots to: plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_subset-refined-modules.pdf\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
