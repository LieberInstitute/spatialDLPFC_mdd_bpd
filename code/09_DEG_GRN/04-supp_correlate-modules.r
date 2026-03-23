setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(dplyr)
  library(pheatmap)
  library(gridExtra)
})
set.seed(123)

# load in GRN adjacency output for correlations
lg.mask = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr.csv")

# load in module sets
#modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_interaction-modules.csv")
modules2 = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules.csv")

# empty input matrix
each_module = unique(modules2$TF)
input.mtx = matrix(NA, nrow=length(each_module), ncol=length(each_module), dimnames = list(each_module, each_module))

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

write.csv(input.mtx, "processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_DEG-modules_target-pearson-correlation.csv")
cat("\nSaved pairwise target correlation matrix for all DEG modules to: processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_DEG-modules_target-pearson-correlation.csv\n")

# plot heatmap and highlight top predictor subsets
mod_subset = c("EIF1",#"COX4I1",
               "UQCRH","PRKAR1A","CAMK2N1","GRIN1","GAD1",
               "MALAT1","A2M","ADAMTS1","IFITM3","CD74",
               "PLP1","HSPA1A","APOE","MT1X")

col_annot = data.frame("is_top"= as.character(each_module %in% mod_subset), row.names=each_module)
annot_colors= list("is_top"=c("FALSE"="white", "TRUE"="black"))
phm = pheatmap(input.mtx, clustering_method="ward.D2",
	breaks= seq(from = -1, to = 1, length.out = 101), 
	treeheight_col = 20, treeheight_row = 20, angle_col=90, fontsize=6,
	annotation_col=col_annot, annotation_row=col_annot, annotation_colors=annot_colors, 
	annotation_legend = FALSE, annotation_names_row = FALSE, annotation_names_col = FALSE,
	main="Pearson corr. of DEG module importance")

pdf(file="plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_DEG-modules_target-pearson-correlation_heatmap.pdf", height=8, width=8)
plot(phm[[4]])
dev.off()

cat("\nHeatmap of all DEG module correlations saved to: plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_DEG-modules_target-pearson-correlation_heatmap.pdf\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()

