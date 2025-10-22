library(dplyr)
library(enrichR)
library(pheatmap)

set.seed(123)

setEnrichrSite("Enrichr") # Human genes
dbs <- listEnrichrDbs()

source("code/08_LR-DEG_analysis/fgsea_functions.r")

gmt_db = "WikiPathways"
res_file = "smoothed-k9-1663"
clust_levels = c("L1","L2","L3.4","L5","L6","WM")
names(clust_levels) = clust_levels

#need full gmt list
if(gmt_db=="Reactome") gmt = .read_gmt("Reactome_2022")
if(gmt_db=="WikiPathways") gmt = .read_gmt("WikiPathways_2024_Human")

gmt.list = group_by(gmt, term) %>% summarise(gene=list(gene)) %>%
  tibble::deframe()

#load fgsea results
fgsea.results = readRDS(paste0("processed-data/08_LR-DEG_analysis/", res_file, "_", gmt_db, "_fgsea-list.rda"))

#load DEG lists
lrList = readRDS("processed-data/08_LR-DEG_analysis/LR-paired_rev-gene-input_logFC-0.3_lists.rds")
saveList <- readRDS("processed-data/07_dx_DE/LA-LR-overlap_rev-gene-input_lists.rds")

#in future will use groupList as input to lapply to loop over all
groupList = names(fgsea.results)
names(groupList) <- groupList
#start with testing 1 group
x = "NTC.BPD_M"

ledge.df = do.call(rbind, lapply(clust_levels, function(y) {
  lr.both.dir = c(lrList[[paste0(x,"_dn")]][["LR_sig"]][[unlist(strsplit(res_file, "-"))[[1]]]][[y]],
                   lrList[[paste0(x,"_up")]][["LR_sig"]][[unlist(strsplit(res_file, "-"))[[1]]]][[y]])
  la.both.dir = c(saveList[["sig_genes"]][[paste0(x,"_dn")]][[paste0("adj_", substr(res_file, start=0, stop=2))]],
                  saveList[["sig_genes"]][[paste0(x,"_up")]][[paste0("adj_", substr(res_file, start=0, stop=2))]])
  
  ledge = tibble::deframe(fgsea.results[[x]][[y]][,c("pathway","leadingEdge")])
  
  d1 = data.frame("term"=names(ledge),
                  "NES"=fgsea.results[[x]][[y]]$NES,
                  "LR_in_leadingEdge"=sapply(lapply(ledge, intersect, y=lr.both.dir), length), 
                  "LA_in_leadingEdge"=sapply(lapply(ledge, intersect, y=la.both.dir), length),
                  "LR.or.LA_in_leadingEdge"=sapply(lapply(ledge, intersect, y=union(la.both.dir, lr.both.dir)), length), 
                  "source"=rep(unlist(strsplit(res_file, "-"))[[1]], length(ledge)),
                  "cluster"=rep(y, length(ledge)),
                  row.names=NULL)
  return(d1)
}))
rownames(ledge.df) <- NULL

## pretty names
if(gmt_db=="Reactome") {
  p.names = sapply(strsplit(ledge.df$term, " R-HSA-"), function(x) x[[1]])
} 
if(gmt_db=="WikiPathways") {
  p.names = sapply(strsplit(ledge.df$term, " WP"), function(x) x[[1]])
}
p.names = sapply(p.names, function(x) {
  c1 = unlist(strwrap(x, width=75))
  if(length(c1)>1) {
    return(paste(c1[[1]],"[...]"))
  } else {
    return(c1[[1]])
  }
})
ledge.df$name = p.names

#pick a subset: should at least keep WM and M.V separate because of greatly increased number of results
tmp1 = filter(ledge.df, cluster %in% c("L2","L3.4","L2.3","L4","L5","L6")) 

md1 = tidyr::pivot_wider(tmp1[,c("name","cluster","NES")], 
                         names_from="cluster", values_from="NES",
                         values_fill=0)
mm1 = as.matrix(md1[,-1])
rownames(mm1) = md1$name

#color limits for NES
color.limits = ceiling(max(abs(mm1)))
if(color.limits-max(abs(mm1))>.5) color.limits = color.limits-.5
seq.len = color.limits*4
pal.len = seq.len-1

phm = pheatmap(mm1, color=colorRampPalette(rev(RColorBrewer::brewer.pal(n=5, "RdBu")))(pal.len),
               breaks= seq(-color.limits, color.limits, length.out=seq.len),
               fontsize_row = 7, treeheight_col = 12, treeheight_row = 12,
               clustering_method="complete", silent=T, angle_col=45,
               main="NES of sig. terms as fill color (0 if not sig.)")
#plot(phm[[4]])

#now color with the number of LR DEGs in leading edge
md2 = tidyr::pivot_wider(tmp1[,c("name","cluster","LR_in_leadingEdge")], 
                         names_from="cluster", values_from="LR_in_leadingEdge",
                         values_fill=0)
mm2 = as.matrix(md2[,-1])
rownames(mm2) = md2$name
col2 = colorRampPalette(RColorBrewer::brewer.pal(n=5, "PuRd"))(max(mm2, na.rm=T)+1)
col2[1] = "white"
phm2 = pheatmap(mm2[phm$tree_row$label[phm$tree_row$order],
                    phm$tree_col$label[phm$tree_col$order]], 
                color=col2,
                fontsize_row = 7, cluster_rows=F, cluster_cols=F, 
                silent=T, angle_col=45,
                main="Fill: # L-R DEGs")
#plot(phm2[[4]])

#now color with the number of LR or LA DEGs in the leading edge
md3 = tidyr::pivot_wider(tmp1[,c("name","cluster","LR.or.LA_in_leadingEdge")], 
                         names_from="cluster", values_from="LR.or.LA_in_leadingEdge",
                         values_fill=0)
mm3 = as.matrix(md3[,-1])
rownames(mm3) = md3$name
col3 = colorRampPalette(RColorBrewer::brewer.pal(n=5, "PuRd"))(max(mm3, na.rm=T)+1)
col3[1] = "white"
phm3 = pheatmap(mm3[phm$tree_row$label[phm$tree_row$order],
                    phm$tree_col$label[phm$tree_col$order]], 
                color=col3,
                fontsize_row = 7, cluster_rows=F, cluster_cols=F,
                silent=T, angle_col=45,
                main="Fill: # L-R or L-A DEGs")
#plot(phm3[[4]])

pdf(file="plots/08_LR-DEG_analysis/hmp_test.pdf", width=7, height=8)
plot(phm[[4]])
plot(phm2[[4]])
plot(phm3[[4]])
dev.off()
