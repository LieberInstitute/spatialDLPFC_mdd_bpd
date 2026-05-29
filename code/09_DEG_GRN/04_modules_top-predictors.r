setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(igraph)
})
set.seed(123)

# load in GRNBOOST2 adj output
lg.mask = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr.csv")

# load DEGs
source("code/09_DEG_GRN/load_DEGs.r")
mbv.degs = unique(do.call(rbind, sigList)$gene_name)

cat("\nNumber of MBv DEGs (F-test adj. p<.05, t-test adj. p<.05):", length(mbv.degs))
cat("\nNumber of MBv DEGs represented in GRN output:", length(intersect(mbv.degs, lg.mask$target)), "\n")

cat("\nPick stricter rho filter for positive regulation (initially rho > 0.03)\n")
tmp = filter(lg.mask, regulation==1)
summary(tmp$rho)
cat("\nNew threshold: rho > 0.2\n")

# select top predictors first, use importance distribution to determine importance cutoff for modules
top.predictors = filter(lg.mask, rho>.2) %>% group_by(target) %>% slice_max(n=1, importance)
#min(top.predictors$importance) #.05
#min(filter(top.predictors, target %in% mbv.degs)$importance) #.06
#plot(ecdf(filter(top.predictors, target %in% mbv.degs)$importance))
#abline(v=1, col="red")
## min importance filter of 1 should work

cat("\nInteraction module criteria: importance > 1, rho > 0.2, at least 19 genes meeting these critera (20 including self)\n")
modules = filter(lg.mask, importance>1, rho>.2) %>% 
  group_by(regulation, TF) %>% add_tally(name="module_size") %>%
  filter(module_size>=19)

cat("\nNumber of interaction modules:", length(unique(modules$TF)))
cat("\nNumber of DEGs represented in interaction modules:", length(intersect(mbv.degs, modules$target)), "\n")

#write.csv(modules, "processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_interaction-modules.csv", row.names=F)
#cat("\nSaved adjacency output filtered to modules to: processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_interaction-modules.csv\n")

# while we're at it, filter to only DEG targets and DEG modules must have at least 10 DEGs
deg.modules = filter(modules, target %in% mbv.degs) %>% 
  group_by(TF) %>% add_tally(name="n_DEGs") %>%
  filter(n_DEGs>=10)
write.csv(deg.modules, "processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules.csv", row.names=F)
cat("\nSaved adjacency output filtered to modules (with at least 10 DEGs) to: processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules.csv\n\n")

# go back to top predictors
cat("\nSelect top positive predictors for each target...\n")
cat("\nNumber of unique top predictors:", length(unique(top.predictors$TF)))

cat("\nNumber of unique top predictors that are interaction modules:", length(unique(intersect(top.predictors$TF, modules$TF))), "\n")

cat("\nFind top predictor networks for MBv DEGs...\n")
e.df = filter(top.predictors, target %in% mbv.degs) %>% 
#              TF %in% modules$TF) %>% 
  select(from=TF, to=target, weight=importance, rho)

n.df = data.frame("node"=union(e.df$from, e.df$to),
                  "is_DEG"=union(e.df$from, e.df$to) %in% mbv.degs,
                  "is_module"= union(e.df$from, e.df$to) %in% modules$TF,
		  "is_DEG.module"= union(e.df$from, e.df$to) %in% deg.modules$TF)

e.reds = colorRampPalette(RColorBrewer::brewer.pal(n=6, "Reds"))(100)
e.df$rho_scaled= round(scales::rescale(e.df$rho, to=c(20,100)),0)


g <- graph_from_data_frame(e.df, directed=T, vertices=n.df)
set.seed(123)
E(g)$weight <- e.df$weight
E(g)$color = e.reds[e.df$rho_scaled]
E(g)$width = scales::rescale(sqrt(e.df$weight), to=c(1,5))

V(g)$size = ifelse(n.df$is_module, 5, 3)
V(g)$size = ifelse(n.df$is_DEG.module, 7, V(g)$size)
V(g)$label.cex = ifelse(n.df$is_module, 1, .5)
V(g)$label.color = "black"
V(g)$color = ifelse(n.df$is_DEG, "grey80", "white")

saveRDS(g, "processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_top-predictor-networks_igraph.rda")
cat("\nSaved top predictor igraph object to: processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_top-predictor-networks_igraph.rda\n")

# break up into subcomponents 
subg = components(g)
#subg$csize

dg <- decompose(g) 
dgl = lapply(dg, function(x) {
  set.seed(123)
  layout= layout_with_kk(x, weights=sqrt(E(x)$weight))
  return(list("igraph"=x, "layout"=layout, "no_mods"=sum(V(x)$is_TF)==0))
})

#plan pdf layout
pdf.layout = data.frame("key"=seq(length(subg$csize)), "csize"=subg$csize, "no_mods"=sapply(dgl, function(x) x$no_mods))
pdf.layout$own_page = subg$csize>30

# filter to subnetworks with >10 DEGs
pdf.layout= filter(pdf.layout, csize>10)

#saveRDS(dgl[pdf.layout$key], "processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_top-predictor-networks.rda")
#cat("\nSaved top predictor network components igraph and layout list to: processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_top-predictor-networks.rda\n")

cat("\nPlot igraph...\n")
#1 per page
list1 = as.list(pdf.layout$key[pdf.layout$own_page])
#3 per page
list2 = list()
rest = filter(pdf.layout, own_page==F) %>% arrange(desc(csize))
while(nrow(rest)>0) {
  add_list = c(rest$key[[1]], tail(rest, 2)$key)
  list2[[(length(list2)+1)]] = add_list
  rest=filter(rest, !key %in% add_list)
}

#setdiff(pdf.layout$key, c(unlist(list1), do.call(c, list2))) #good, didn't miss any

layout(mat=as.matrix(c(1)), widths=c(1), heights=c(1))
pdf(file="plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_top-predictor-networks.pdf", 
    width=8, height=8)
for(i in list1) {
  plot(dgl[[i]]$igraph, layout=dgl[[i]]$layout, 
       vertex.label.family="sans", vertex.frame.width=0, vertex.frame.color=NA)
}
layout.matrix <- matrix(c(1, 1, 2, 3), nrow = 2, byrow=T)
layout(mat = layout.matrix,
       heights = c(1, 1), # Heights of the two rows
       widths = c(1, 1)) # Widths of the two columns
for(i in list2) {
  for(j in i) {
    plot(dgl[[j]]$igraph, layout=dgl[[j]]$layout, 
         vertex.label.family="sans", vertex.frame.width=0, vertex.frame.color=NA)
  }
}
dev.off()

cat("\nTop predictor PDF saved to: plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_top-predictor-networks.pdf\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
