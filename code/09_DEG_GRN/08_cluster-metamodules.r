setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(pheatmap)
})
set.seed(123)

GENESET = "13844"
MODULES_TYPE = "modules-top20-Fadjp05"
ADJ_TYPE = "top20percent-Fadjp05"


adj <- read.csv(paste0("processed-data/09_DEG_GRN/spe-n119_", GENESET, "-no-lowUMI_adj_with-logcounts-corr_",ADJ_TYPE,".csv"))
aucell_corr = read.csv(paste0("processed-data/09_DEG_GRN/spe-n119_",GENESET,"-no-lowUMI_logcounts_",MODULES_TYPE,"_AUCell-correlation.csv"))

#format aucell
cnames = gsub("\\(\\+\\)", "_act", gsub("\\(\\-\\)", "_rep", aucell_corr[,1]))
aucell_mtx = as.matrix(aucell_corr[,-1])
#grep("rep", colnames(aucell_mtx), value=T) #only ATP1B1_rep
colnames(aucell_mtx) = cnames
rownames(aucell_mtx) = cnames

# initial clustering
phm = pheatmap(aucell_mtx, clustering_method = "ward.D2", fontsize=6, silent=T)

initial_annot = data.frame("h90"=paste0("h",cutree(phm$tree_col, h=quantile(phm$tree_col$height, 
                                                                        probs=c(.90)))),
                       row.names=phm$tree_col$labels)

df1 = data.frame("clust"="initial",
                 "hcut"=.9,
                 "hID"=c("h11","h12","h8","h4","h7","h9","h6","h3"),
                 "mmID"=c("mm1","mm2","mm3","mm4","mm5","mm6","mm7","mm8"),
                 "mmDesc"=c("myelin","translation","housekeeping","eif1/aplp1","mitochondrial matrix?","electron transport","neuron housekeeping?","neuronal function"))
df1$modules = sapply(df1$hID, function(x) {
  m1 = rownames(initial_annot)[initial_annot$h90==x] 
  paste(gsub("_act|_rep","",m1), collapse="/")
})

# split h3 to separate KIF5C_act and ATP1B1_rep
df1[df1$mmID=="mm8","modules"] = paste(setdiff(unlist(strsplit(df1$modules[df1$mmID=="mm8"], "/")), c("KIF5C","ATP1B1")), collapse="/")
df1 = mutate(df1, dir=1) %>% add_row(clust="initial", hcut=.9, hID="h3", mmID=c("mm9","mm10"), mmDesc=c("KIF5C_act","ATP1B1_rep"), modules=c("KIF5C","ATP1B1"), dir=c(1,-1))


phm2 = pheatmap(aucell_mtx[rownames(initial_annot)[!initial_annot$h90 %in% c("h11","h12","h8","h4","h9","h7","h6","h3")],],
                clustering_method="ward.D2", silent=T)
                #annotation_row = col_annot, annotation_col= col_annot, annotation_colors = annot_colors, fontsize_col=6,
                #main="rest", 

sub_annot = data.frame("h75"=paste0("h",cutree(phm2$tree_row, h=quantile(phm2$tree_row$height, 
                                                                         probs=c(.75)))),
                       row.names=phm2$tree_row$labels)

df2 = data.frame("clust"="sub",
                 "hcut"=.75,
                 "hID"=c("h6","h4","h8","h1","h2","h9","h10","h3","h11",
                         "h7","h5","h13","h12"),
                 "mmID"=c("mm11","mm12","mm13","mm14","mm15","mm16","mm17","mm18","mm19",
                          "mm20","mm21","mm22","mm23"),
                 "mmDesc"=c("microglia","astrocyte","GABAergic","vascular","ECM?","glial activation?","heatshock","BBB?","batch effect",
                            "FOS","AQP1","SERPINA3","SCD5"))
df2$modules = sapply(df2$hID, function(x) {
  m1 = rownames(sub_annot)[sub_annot$h75==x]
  paste(gsub("_act|_rep","",m1), collapse="/")
})
df2$dir = 1

mm.df = rbind(df1, df2)
write.csv(mm.df, paste0("processed-data/09_DEG_GRN/spe-n119_",GENESET,"-no-lowUMI_", MODULES_TYPE, "_to-metamodules_key.csv"), row.names=F)
cat("\nSaved modules-to-metamodules key table to:", paste0("processed-data/09_DEG_GRN/spe-n119_",GENESET,"-no-lowUMI_", MODULES_TYPE, "_to-metamodules_key.csv"),"\n")

mm.long = tidyr::separate_rows(mm.df, modules, sep="/") %>%
  mutate(modules_key=paste0(modules, as.character(factor(dir, levels=c(-1,1), labels=c("_rep","_act")))))


mod_annot = data.frame("metaMod"=mm.long$mmID, row.names=mm.long$modules_key)
cp = c("#821010","#9E8810","#2C822C","#21609E","#826168",
       "red3","goldenrod1","limegreen","dodgerblue","pink3",
       "pink", "#E9D878", "#8ED38E", "#85B8E9",
       "#A1591D","black","#D28015","#6B7B2E","#615D94",
       "#F3DFC4", "grey","#E8BF8A", "turquoise4")
names(cp) = names(table(mod_annot[,1]))
annot_colors=list("metaMod"=cp)

pdf(file = paste0("plots/09_DEG_GRN/spe-n119_",GENESET,"-no-lowUMI_", MODULES_TYPE, "_to-metamodules_heatmap.pdf"))
pheatmap(aucell_mtx, clustering_method = "ward.D2", fontsize=6, 
         treeheight_row = 12, treeheight_col = 20,
         show_rownames=F, show_colnames=F, annotation_names_row = F, annotation_names_col = F, 
         annotation_col = mod_annot, annotation_row = mod_annot, annotation_colors = annot_colors)
dev.off()
cat("\nHeatmap of metamodules saved to:", paste0("plots/09_DEG_GRN/spe-n119_",GENESET,"-no-lowUMI_", MODULES_TYPE, "_to-metamodules_heatmap.pdf"), "\n")


#metamodules to targets
mm.adj = do.call(rbind, lapply(mm.df$mmID, function(x) {
  mm.filt = filter(mm.long, mmID==x)
  adj.filt = filter(adj, TF %in% mm.filt$modules, regulation==unique(mm.filt$dir))
  
  missingTF= setdiff(mm.filt$modules, adj.filt$target)
  if(length(missingTF)>0) {
    if(nrow(mm.filt)>1) {
      warning(paste0(x,": Some TFs missing from targets! (not singlular and therefore not added)"))
    } else {
      if(mm.filt$dir<0) {
        warning(paste0(x,": TF missing from repressing module, not added"))
      } else {
        adj.filt = add_row(adj.filt, TF="any", target=missingTF, importance=1, regulation=1, rho=1)
      }
    }
  }
  adj.summ = group_by(adj.filt, target, regulation) %>% 
    summarise(n_modules=n(), importance=mean(importance), rho=mean(rho)) %>%
    mutate(TF=x, n_TFs=nrow(mm.filt))
  return(as.data.frame(adj.summ)[,c("TF","target","importance","regulation","rho","n_TFs","n_modules")])
}))
#mm15 RHOB is missing and RHOB is well predicted by others so not adding


##consider filtering to genes that were predicted in 30% of modules?
#head(mm.adj)
#print(n=23, group_by(mm.adj, TF) %>% tally() %>% arrange(TF))
#print(n=23, mutate(mm.adj, prop_TF_modules=n_modules/n_TFs) %>% filter(prop_TF_modules>.3) %>%
#  group_by(TF) %>% tally())

mm.adj$n_TFs <- NULL
mm.adj$n_modules <- NULL

write.csv(mm.adj, paste0("processed-data/09_DEG_GRN/spe-n119_",GENESET,"-no-lowUMI_adj_with-logcounts-corr_", ADJ_TYPE, "-metaModules.csv"), row.names=F)
cat("\nMetamodule adjacency file saved to:", paste0("processed-data/09_DEG_GRN/spe-n119_",GENESET,"-no-lowUMI_adj_with-logcounts-corr_", ADJ_TYPE, "-metaModules.csv"), "\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
