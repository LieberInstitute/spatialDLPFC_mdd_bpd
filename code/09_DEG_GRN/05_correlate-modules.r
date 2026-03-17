setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(dplyr)
  library(pheatmap)
  library(gridExtra)
})
set.seed(123)

GENESET = "13162"
MODULES_TYPE = "modules-top20-Fadjp05"
ADJ_TYPE = "top20percent-Fadjp05"
cat("\nGene set:", GENESET, "\n")
cat("\nAdjacency type:", ADJ_TYPE, "\n")

adj_raw <- read.csv(paste0("processed-data/09_DEG_GRN/spe-n119_", GENESET, "-no-lowUMI_adj_with-logcounts-corr.csv")) %>%
  mutate(modules_key=paste0(TF, as.character(factor(regulation, levels=c(-1,1), labels=c("_rep","_act")))))
adj <- read.csv(paste0("processed-data/09_DEG_GRN/spe-n119_", GENESET, "-no-lowUMI_adj_with-logcounts-corr_",ADJ_TYPE,".csv")) %>%
  mutate(modules_key=paste0(TF, as.character(factor(regulation, levels=c(-1,1), labels=c("_rep","_act")))))


modules = unique(adj$modules_key)
cat("\nNumber of modules:", length(modules), "\n")

input.mtx = matrix(NA, nrow=length(modules), ncol=length(modules), dimnames = list(modules, modules))
cat("\nPearson correlation between module targets...\n")
for(i in 1:(nrow(input.mtx)-1)) {
  complete.corr = colnames(input.mtx)[(i+1):ncol(input.mtx)]
  i=rownames(input.mtx)[i]
  for(j in complete.corr) {
    #extract target genes and exclude selves 
    tmp1 = filter(adj, modules_key==i)#,
    #              target %in% mbv.degs2) # ADDED THIS FILTER FOR DEG MTX
    i1 = unique(tmp1$TF)
    tmp2 = filter(adj, modules_key==j)#,
    #              target %in% mbv.degs2) # ADDED THIS FILTER FOR DEG MTX
    j1 = unique(tmp2$TF)
    geneList = setdiff(union(tmp1$target, tmp2$target),
                       c(i1,j1))
    #pull the raw adj (before any filtering) and format for correlation
    check = filter(adj_raw, modules_key %in% c(i,j), target %in% geneList) %>%
      mutate(modules_key=factor(modules_key, levels=c(i,j), labels=c("I","J"))) %>%
      select(modules_key, target, importance) %>% tidyr::pivot_wider(names_from="modules_key", values_from="importance", values_fill=0)
    corrij= cor(check$I, check$J)
    input.mtx[i,j] = corrij
    input.mtx[j,i] = corrij
  }
}
write.csv(input.mtx, paste0("processed-data/09_DEG_GRN/spe-n119_",GENESET,"-no-lowUMI_logcounts_",MODULES_TYPE,"_target-pearson-correlation.csv"))
cat("\nCorrelation matrix saved to:", paste0("processed-data/09_DEG_GRN/spe-n119_",GENESET,"-no-lowUMI_logcounts_",MODULES_TYPE,"_target-pearson-correlation.csv"),"\n")

phm = pheatmap(input.mtx, clustering_method = "ward.D2", silent=T)
initial_annot = data.frame("h95"=paste0("h",cutree(phm$tree_col, h=quantile(phm$tree_col$height,
                                                                            probs=c(.95)))),
                           row.names=phm$tree_col$labels)
cat("\nModule target correlation clusters (at 95% height):\n")
table(initial_annot[,1])
write.csv(initial_annot, paste0("processed-data/09_DEG_GRN/spe-n119_",GENESET,"-no-lowUMI_logcounts_",MODULES_TYPE,"_target-pearson-correlation_annotation.csv"))

cp = c("red3","goldenrod1","limegreen","dodgerblue","pink3",
       "#821010","#9E8810","#2C822C","#21609E","#826168",
       "pink", "#E9D878", "#8ED38E", "#85B8E9",
       "#A1591D","black","#D28015","#6B7B2E","#615D94",
       "#F3DFC4", "grey","#E8BF8A", "turquoise4")
cp2 = cp[1:length(unique(initial_annot[,1]))]
names(cp2) = names(table(initial_annot[,1]))
annot_colors=list("h95"=cp2)


phm1 = pheatmap(input.mtx, clustering_method = "ward.D2", fontsize=8,
	 breaks= seq(from = -1, to = 1, length.out = 101), silent=T,
         treeheight_row = 12, treeheight_col = 20, border_color=NA,
         show_rownames=F, show_colnames=F, annotation_names_row = F, annotation_names_col = F,
         annotation_col = initial_annot, annotation_row = initial_annot, annotation_colors = annot_colors)

corder = phm1$tree_col$labels[phm1$tree_col$order]
clist = unique(initial_annot[corder,1])

plist <- lapply(clist, function(x){
  metamod = intersect(corder, rownames(initial_annot)[initial_annot[,1]==x])
  pp = pheatmap(input.mtx[metamod, metamod], cluster_rows=F, cluster_cols=F,
                breaks= seq(from = -1, to = 1, length.out = 101),
                fontsize=8, angle_col=90, main=x, #silent=T,
                annotation_names_row = F, annotation_names_col = F,
                annotation_col = initial_annot, annotation_row = initial_annot, annotation_colors = annot_colors)
  return(pp[[4]])
})


pdf(file=paste0("plots/09_DEG_GRN/spe-n119_",GENESET,"-no-lowUMI_logcounts_",MODULES_TYPE,"_target-pearson-correlation_heatmap.pdf"))
plot(phm1[[4]])
marrangeGrob(plist, ncol=1, nrow=1, top=NULL)
dev.off()

cat("\nHeatmap saved to:", paste0("plots/09_DEG_GRN/spe-n119_",GENESET,"-no-lowUMI_logcounts_",MODULES_TYPE,"_target-pearson-correlation_heatmap.pdf"),"\n")
cat("Heatmap annotation key saved to:",paste0("processed-data/09_DEG_GRN/spe-n119_",GENESET,"-no-lowUMI_logcounts_",MODULES_TYPE,"_target-pearson-correlation_annotation.csv"),"\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()

