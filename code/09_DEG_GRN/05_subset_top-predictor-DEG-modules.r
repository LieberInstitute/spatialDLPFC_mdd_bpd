setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(dplyr)
	library(pheatmap)
	library(ggplot2)
	library(igraph)
	library(SpatialExperiment)
})

set.seed(123)

# load in GRN adjacency output for correlations
lg.mask = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr.csv")

# load in module sets
modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_interaction-modules.csv")
modules2 = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules.csv")

#pick out top predictor DEG network subset
mod_subset = c("COX4I1","UQCRH","PRKAR1A","CAMK2N1","FAIM2","GAD1",
               "SNHG14","APOE","A2M","ADAMTS1","IFITM3","CD74",
               "FTL","APLP1","PLP1","HSPA1A","MT1X")
cat("\nNumber of INITIAL top predictor DEG network modules to highlight:", length(mod_subset), "\n")


input.mtx = matrix(NA, nrow=length(mod_subset), ncol=length(mod_subset), dimnames = list(mod_subset, mod_subset))
for(i in 1:(nrow(input.mtx)-1)) {
  complete.corr = colnames(input.mtx)[(i+1):ncol(input.mtx)]
  i=rownames(input.mtx)[i]
  for(j in complete.corr) {
    tmp1 = filter(modules, TF==i)
    tmp2 = filter(modules, TF==j)
    geneList = setdiff(union(tmp1$target, tmp2$target),
                       c(i,j))
    #pull the raw adj (before any filtering) and format for correlation
    check = filter(lg.mask, TF %in% c(i,j), target %in% geneList) %>%
      mutate(modules_key=factor(TF, levels=c(i,j), labels=c("I","J"))) %>%
      select(modules_key, target, importance) %>% tidyr::pivot_wider(names_from="modules_key", values_from="importance", values_fill=0)
    corrij= cor(check$I, check$J)
    input.mtx[i,j] = corrij
    input.mtx[j,i] = corrij
  }
}

## identify highly correlated DEG modules and pick one
#input.mtx[matrixStats::rowMaxs(input.mtx, na.rm=T)>.2,
#          matrixStats::colMaxs(input.mtx, na.rm=T)>.2]
# keep UQCRH
cat("\nTwo top predictor networks represent highly correlated target gene modules...")
cat("\nRemove COX4I1 and keep UQCRH...\n")
exclude = c("COX4I1")
mod_subset = setdiff(mod_subset, exclude)

cat("\nNumber of FINAL top predictor DEG network modules to highlight:", length(mod_subset), "\n")
cat("\nTop predictor DEG network modules to highlight:\n")
print(mod_subset)


# filter DEG modules to top predictor subset
tmp = filter(modules2, TF %in% mod_subset)

source("code/09_DEG_GRN/load_DEGs.r")
mbv.degs = unique(do.call(rbind, sigList)$gene_name)

cat("\nNumber of DEGs present in GRN output:", length(intersect(mbv.degs, lg.mask$target)),"\n")
cat("\nNumber of DEGs present in the", length(unique(modules$TF)), "interaction modules:", length(intersect(mbv.degs, modules$target)),"\n")
cat("\nNumber of DEGs present in the", length(unique(modules2$TF)), "DEG modules:", length(intersect(mbv.degs, modules2$target)),"\n")
cat("\nNumber of DEGs present in the", length(mod_subset), "top predictor DEG modules:", length(intersect(mbv.degs, tmp$target)),"\n\n")

# plot top predictor modules heatmaps
phm1 = pheatmap(input.mtx[mod_subset,mod_subset], clustering_method="ward.D2",
                   breaks= seq(from = -1, to = 1, length.out = 101), 
                   treeheight_col = 20, treeheight_row = 20, angle_col=90,
                main="Pearson corr. of module importance (all targets)")

# second heatmap based on correlations of only the DEGs in the modules
input.mtx = matrix(NA, nrow=length(mod_subset), ncol=length(mod_subset), dimnames = list(mod_subset, mod_subset))
for(i in 1:(nrow(input.mtx)-1)) {
  complete.corr = colnames(input.mtx)[(i+1):ncol(input.mtx)]
  i=rownames(input.mtx)[i]
  for(j in complete.corr) {
    tmp1 = filter(modules2, TF==i)
    tmp2 = filter(modules2, TF==j)
    geneList = setdiff(union(tmp1$target, tmp2$target),
                       c(i,j))
    #pull the raw adj (before any filtering) and format for correlation
    check = filter(lg.mask, TF %in% c(i,j), target %in% geneList) %>%
      mutate(modules_key=factor(TF, levels=c(i,j), labels=c("I","J"))) %>%
      select(modules_key, target, importance) %>% tidyr::pivot_wider(names_from="modules_key", values_from="importance", values_fill=0)
    corrij= cor(check$I, check$J)
    input.mtx[i,j] = corrij
    input.mtx[j,i] = corrij
  }
}

phm2 = pheatmap(input.mtx, clustering_method="ward.D2",
                   breaks= seq(from = -1, to = 1, length.out = 101), 
                   treeheight_col = 20, treeheight_row = 20, angle_col=90,
                main="Pearson corr. of module importance (DEG targets only)")


# plot heatmap of expression for top predictor DEG modules
#load in sce for heatmap and dotplots
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-no-lowUMI-dotplot_dx-sex-seurat-pc30.Rdata")
cond_sex = c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M")
spe_summ$seurat_label = factor(spe_summ$seurat_qual.genes_pc30.kweight50, levels=c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"),
                               labels=c("M.V","Ast","L2.3","L4","Inb","L5","L6","Olg"))
seurat_levels= c("M.V","Ast","L2.3","L4","Inb","L5","L6","Olg")
spe_summ$sample_id = factor(paste(spe_summ$condition, spe_summ$sex, spe_summ$seurat_label),
                            levels=as.character(outer(cond_sex, seurat_levels, paste)))
colnames(spe_summ) <- spe_summ$sample_id

#load in dotplotDF function and format dataframe for dotplot
source("code/06_pseudobulk/custom_functions.r")

sm.df = dotplotDF(spe_summ, mod_subset, swap_rownames="gene_name",
                  summarize_groups=F, row_data=NULL) %>%
  mutate(gene_name_f=factor(gene_name, levels=rev(mod_subset)),
         clusters=factor(clusters, levels=levels(spe_summ$sample_id), labels=gsub(" ","\n", levels(spe_summ$sample_id))))
sm.df$condition = factor(substr(sm.df$clusters, start=0, stop=3), levels=c("NTC","MDD","BPD")) 
sm.df$sex = factor(substr(sm.df$clusters, start=5, stop=5), levels=c("F","M"))

#make labels
xlbs = levels(sm.df$clusters)
## keep only the odd numbered dx indicators
odd1 = seq(1, length(xlbs), by=2)
xlbs[-odd1] = substr(xlbs[-odd1], start=4, stop=11)

p1 <- ggplot(sm.df, aes(x=clusters, y=gene_name_f))+
  geom_count(aes(shape=sex, color=mean_expr_scaled, size=prop_spots))+
  scale_shape_manual(values=c(20,18))+
  scale_color_gradient(low="white", high="black")+
  scale_x_discrete(labels=xlbs)+
  scale_size(range=c(1,5), limits=c(0,1), breaks=c(0,.5,1))+
  guides(size = guide_legend(override.aes = list(shape = 20)),
         shape = guide_legend(overrisde.aes = list(size=5)))+
  labs(color="Avg. expr.\n(logcounts,\nscaled)", size="Prop. of\nspots",
       y="top predictors", title="Top predictor DEG modules")+
  theme_minimal()+theme(axis.title.x=element_blank(), axis.text.x=element_text(size=7, hjust=0),
                        axis.text.y=element_text(face="italic"), legend.key.size=unit(15,"pt"))




# plot igraph of top predictor DEG modules to highlight

plotModuleNetwork <- function(source_DF, filter_set="none", edge_color="grey", vertex_color="lightgoldenrod1") {
  e.df = select(ungroup(tmp), from=TF, to=target, weight=importance)
  if(length(filter_set)>1) {
    e.df = filter(e.df, from %in% filter_set)
  }
  n.df = data.frame("node"=union(e.df$from, e.df$to),
                    "is_TF"= union(e.df$from, e.df$to) %in% e.df$from)
  
  g <- graph_from_data_frame(e.df, directed=T, vertices=n.df)
  E(g)$weight = e.df$weight
  E(g)$arrow.size = 0
  E(g)$color = edge_color
  
  V(g)$size = ifelse(n.df$is_TF, 5, 1)
  V(g)$label = ifelse(n.df$is_TF, n.df$node, "")
  V(g)$label.color = "black"
  V(g)$color = vertex_color
  
  V(g)$label.family="sans"
  V(g)$frame.width=0
  V(g)$frame.color=NA
  
  set.seed(123)
  lay1= layout_with_fr(g, weights=E(g)$weight)
  
  return(list("igraph"=g, "layout"=lay1))
}

## all
filterSets = list("none", 
                  c("UQCRH","CAMK2N1","PRKAR1A","FAIM2","GAD1"),
                  c("APOE","A2M","ADAMTS1","IFITM3","MT1X","CD74"),
                  c("APLP1","FTL","PLP1","SNHG14","HSPA1A"))
names(filterSets) <- c("all","neuron","bbb","other")

plotList <- lapply(filterSets, function(x) plotModuleNetwork(tmp, filter_set=x))

saveRDS(plotList, "processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_top-predictor-DEG-modules.rda")
cat("\nSaved list of igraphs and layouts for top predictor DEG modules to: processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_top-predictor-DEG-modules.rda\n")

pdf(file="plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_top-predictor-DEG-modules.pdf", width=8, height=8)
layout.matrix <- matrix(c(1, 0, 0, 2), nrow = 2, byrow=T)
layout(mat = layout.matrix,
       heights = c(1, 1), # Heights of the two rows
       widths = c(1, 1)) # Widths of the two columns
plot(phm1[[4]])
plot(phm2[[4]])
layout.matrix <- matrix(c(1), nrow = 1, ncol=1)
layout(mat = layout.matrix,
       heights = c(1), # Heights of the two rows
       widths = c(1)) # Widths of the two columns
p1 #dotplot of expression
for(i in 1:length(plotList)) {
  plot(plotList[[i]]$igraph, layout=plotList[[i]]$layout, main=names(plotList)[[i]])
}
dev.off()
cat("\nSaved all plots to: plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_top-predictor-DEG-modules.pdf\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
