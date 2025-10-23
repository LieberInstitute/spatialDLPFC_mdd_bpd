args = commandArgs(TRUE)
x= args[[1]]
print(x)

setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(enrichR)
	library(pheatmap)
	library(igraph)
	library(gridExtra)
})

set.seed(123)

setEnrichrSite("Enrichr") # Human genes
dbs <- listEnrichrDbs()

source("code/08_LR-DEG_analysis/fgsea_functions.r")
source("code/08_LR-DEG_analysis/jc-igraph_functions.r")

cpList = readRDS("plots/colorPalettes.rds")

#load DEG lists
lrList = readRDS("processed-data/08_LR-DEG_analysis/LR-paired_rev-gene-input_logFC-0.3_lists.rds")
saveList <- readRDS("processed-data/07_dx_DE/LA-LR-overlap_rev-gene-input_lists.rds")

#gmt_db = "Reactome"
gmt_db = "WikiPathways"

#(res_file = "smoothed-k9-1663")
#(clust_levels = c("L1","L2","L3.4","L5","L6","WM"))
#clust_subset = c("L2","L3.4","L5","L6")

(res_file = "seurat-pc30")
(clust_levels = c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))
clust_subset = c("L2.3","L4","Inhb","L5","L6")

names(clust_levels) = clust_levels


#load fgsea results
fgsea.results = readRDS(paste0("processed-data/08_LR-DEG_analysis/", res_file, "_", gmt_db, "_fgsea-list.rda"))

ledge.df = do.call(rbind, lapply(clust_levels, function(y) {
  #pull leading edge genes
  ledge = tibble::deframe(fgsea.results[[x]][[y]][,c("pathway","leadingEdge")])
  
  #pull LR genes
  lr.both.dir = c(lrList[[paste0(x,"_dn")]][["LR_sig"]][[unlist(strsplit(res_file, "-"))[[1]]]][[y]],
                  lrList[[paste0(x,"_up")]][["LR_sig"]][[unlist(strsplit(res_file, "-"))[[1]]]][[y]])
  
  #pull LA genes
  la.both.dir = c(saveList[["sig_genes"]][[paste0(x,"_dn")]][[paste0("adj_", substr(res_file, start=0, stop=2))]],
                  saveList[["sig_genes"]][[paste0(x,"_up")]][[paste0("adj_", substr(res_file, start=0, stop=2))]])
  
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

# heatmaps

# focus on neuron dominant groups
tmp1 = filter(ledge.df, cluster %in% clust_subset,
              LR.or.LA_in_leadingEdge>0) 

#WikiPathways has some identically named terms when you drop the WP ID so need to fix that or else throws error
check1 = distinct(tmp1, term, name) %>% group_by(name) %>% add_tally() %>% filter(n>1)
if(nrow(check1)>0) {
  if(gmt_db=="Reactome") {
    p.names2 = paste0("R-HSA-", sapply(strsplit(check1$term, " R-HSA-"), function(x) x[[2]]))
  } 
  if(gmt_db=="WikiPathways") {
    p.names2 = paste0("WP", sapply(strsplit(check1$term, " WP"), function(x) x[[2]]))
    p.names2 = paste(p.names2, check1$name)
  }
  p.names2 = sapply(p.names2, function(x) {
    c1 = unlist(strwrap(x, width=75))
    if(length(c1)>1) {
      return(paste(c1[[1]],"[...]"))
    } else {
      return(c1[[1]])
    }
  })
  check1$name = p.names2
  #change pretty name to make distinct
  for(i in 1:nrow(check1)) {
    tmp1[tmp1$term==check1[["term"]][i],"name"] = check1[["name"]][i]
  }
}

# heatmap of NES
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
               #fontsize_row = 7, 
               clustering_method="complete", silent=T, angle_col=0,
               cluster_cols=F, treeheight_row = 12, #treeheight_col = 12,
               main="NES of sig. terms as fill color (0 if not sig.)")


# heatmap of the number of LR or LA DEGs in the leading edge
md2 = tidyr::pivot_wider(tmp1[,c("name","cluster","LR.or.LA_in_leadingEdge")], 
                         names_from="cluster", values_from="LR.or.LA_in_leadingEdge",
                         values_fill=0)
mm2 = as.matrix(md2[,-1])
rownames(mm2) = md2$name
col2 = colorRampPalette(RColorBrewer::brewer.pal(n=5, "PuRd"))(max(mm2, na.rm=T)+1)
col2[1] = "white"
phm2 = pheatmap(mm2[phm$tree_row$label[phm$tree_row$order],],
                    #phm$tree_col$label[phm$tree_col$order]], 
                color=col2, cluster_rows=F, cluster_cols=F,
                #fontsize_row = 7, 
                silent=T, angle_col=0,
                main="Fill: # L-R or L-A DEGs")


# now switch to looking at genes
ledge.df2 = do.call(rbind, lapply(clust_levels, function(y) {
  #expand df by leading edge genes
  ledge = select(fgsea.results[[x]][[y]], term=pathway, NES, gene_name=leadingEdge2) %>% 
    tidyr::separate_rows(gene_name, sep="/")
  
  #pull LR genes
  lr.both.dir = c(lrList[[paste0(x,"_dn")]][["LR_sig"]][[unlist(strsplit(res_file, "-"))[[1]]]][[y]],
                  lrList[[paste0(x,"_up")]][["LR_sig"]][[unlist(strsplit(res_file, "-"))[[1]]]][[y]])
  ledge = mutate(ledge, is_LR = gene_name %in% lr.both.dir)
  #pull LA genes
  la.both.dir = c(saveList[["sig_genes"]][[paste0(x,"_dn")]][[paste0("adj_", substr(res_file, start=0, stop=2))]],
                  saveList[["sig_genes"]][[paste0(x,"_up")]][[paste0("adj_", substr(res_file, start=0, stop=2))]])
  ledge = mutate(ledge, is_LA = gene_name %in% la.both.dir)
  
  mutate(ledge, is_LR.or.LA= is_LR|is_LA, cluster=y)
}))  

tmp2 = left_join(tmp1, ledge.df2, by=c("term","cluster","NES")) %>%
  filter(is_LR.or.LA==T)

write.csv(tmp2, file=paste0("processed-data/08_LR-DEG_analysis/", gsub("\\.","-", gsub("_","-",x)), "_", res_file, "_", gmt_db, "_neuronal-fgsea_filtered-DEG-leadingeEdge.csv"), row.names=F)
cat("\nFiltered GSEA data frame (expanded by genes) saved to:", paste0("processed-data/08_LR-DEG_analysis/", gsub("\\.","-", gsub("_","-",x)), "_", res_file, "_", gmt_db, "_neuronal-fgsea_filtered-DEG-leadingeEdge.csv"), "\n")


tmp3 = group_by(tmp2, term, name, gene_name) %>% tally(name="n_clusters") %>%
  tidyr::pivot_wider(names_from="gene_name", values_from="n_clusters", values_fill=0)
m3 = as.matrix(tmp3[,3:ncol(tmp3)])
rownames(m3) <- tmp3$name
col3 = colorRampPalette(RColorBrewer::brewer.pal(n=5, "YlOrBr"))(max(m3, na.rm=T)+1)
col3[1] = "white"
phm3 = pheatmap(m3[phm$tree_row$label[phm$tree_row$order],],
         color=col3, treeheight_col = 12, silent=T, 
         cluster_rows=F, fontsize_row=7, fontsize_col=7, angle_col=90,
         main="Fill: # clusters")


# igraph 
out2 = formatJaccardIGRAPH(fgsea.results, .group=x, .cluster=clust_subset)
#filter out sig terms without any DEG
out2$nodes = filter(out2$nodes, LR.or.LA_in_leadingEdge>0)
out2$edges = filter(out2$edges, reference %in% out2$nodes$term, query %in% out2$nodes$term)

saveRDS(out2, file=paste0("processed-data/08_LR-DEG_analysis/", gsub("\\.","-", gsub("_","-",x)), "_", res_file, "_", gmt_db, "_neuronal-fgsea_jc-igraph.rda"))
cat("\nNode and edge data frames saved in list to:", paste0("processed-data/08_LR-DEG_analysis/", gsub("\\.","-", gsub("_","-",x)), "_", res_file, "_", gmt_db, "_neuronal-fgsea_jc-igraph.rda"), "\n")

#quick fix:
cpList$seurat.bright = cpList$transfer.bright
outList2 = generateIGRAPH(out2)


#save plots
pdf(file=paste0("plots/08_LR-DEG_analysis/", gsub("\\.","-", gsub("_","-",x)), "_", res_file, "_", gmt_db, "_neuronal-fgsea-summary.pdf"), height=12, width=12)
grid.arrange(phm[[4]], top=paste(gsub("_"," ",x), gmt_db))
grid.arrange(phm2[[4]], top=paste(gsub("_"," ",x), gmt_db))
grid.arrange(phm3[[4]], top=paste(gsub("_"," ",x), gmt_db))
plot(simplify(outList2[["igraph"]]), layout=outList2[["layout"]], 
     edge.width=E(outList2[["igraph"]])$jc*5, 
     vertex.label.family="sans", #vertex.label.font=2, 
     vertex.frame.width=3,
     main=paste(paste0(gsub("_", " ", x)," ", gmt_db, ":"), outList2[["title_text"]]),
     sub=outList2[["sub_text"]]
)
legend("bottomleft", legend=c("Depleted","Depleted (with LR DEG)",
                              "Enriched","Enriched (with LR DEG)"),
       pch=16, pt.cex=1, cex=.7,
       col=c("#CFEBF7","skyblue","#FFC0B5","tomato"))
dev.off()
cat("\nPlots of neuronal GSEA summary saved to:", paste0("plots/08_LR-DEG_analysis/", gsub("\\.","-", gsub("_","-",x)), "_", res_file, "_", gmt_db, "_neuronal-fgsea-summary.pdf"), "\n")


cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
