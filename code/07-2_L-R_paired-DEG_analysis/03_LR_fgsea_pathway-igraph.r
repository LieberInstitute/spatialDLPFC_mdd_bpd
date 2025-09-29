setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(pheatmap)
	library(gridExtra)
	library(grid)
	library(igraph)
})

set.seed(123)

source("code/07-2_L-R_paired-DEG_analysis/fgsea_functions.r")
source("code/07-2_L-R_paired-DEG_analysis/jc-igraph_functions.r")

# load resources
cpList = readRDS("plots/colorPalettes.rds")

coefList = c("NTC.MDD_F",#"NTC.MDD_M",
             "NTC.BPD_F","NTC.BPD_M","MDD.BPD_F","MDD.BPD_M")
names(coefList) = coefList

# make safe legend position list (based on plotting igraph with seed 123)
lpList = list()
lpList[["NTC.MDD_F"]] = c("sm_L2"="bottomleft","sm_L3.4"="bottomleft","se_L2.3"="bottomleft","se_L4"="topleft","sm_L5"="bottomleft","se_L5"="topleft","sm_L6"="bottomleft","se_L6"="topleft")
lpList[["NTC.BPD_F"]] = c("sm_L2"="bottomleft","sm_L3.4"="bottomleft","se_L2.3"="bottomleft","se_L4"="topleft","sm_L5"="bottomleft","se_L5"="bottomleft","sm_L6"="topleft","se_L6"="bottomleft")
lpList[["NTC.BPD_M"]] = c("sm_L2"="bottomleft","sm_L3.4"="left","se_L2.3"="bottomleft","sm_L5"="topleft","se_L5"="bottomleft","sm_L6"="topleft","se_L6"="bottomleft")
lpList[["MDD.BPD_F"]] = c("sm_L2"="bottomright","se_L2.3"="topleft","se_L6"="topleft")
lpList[["MDD.BPD_M"]] = c("sm_L2"="bottomleft","sm_L3.4"="bottomleft","se_L2.3"="bottomleft","se_L4"="bottomleft","sm_L5"="bottomleft","se_L5"="left","se_L6"="bottomleft")

# set gsea results to plot
gmt_db = "Reactome"
#gmt_db = "WikiPathways"

for(z in coefList) {
  cat("\n\n************\n")
  cat(z,"\n")
  
  gmtList <- readRDS(paste0("processed-data/07-2_L-R_paired-DEG_analysis/", z, "_", gmt_db, "_GSEA.rds"))
  cat("\nLoad saved reactome results from:", paste0("processed-data/07-2_L-R_paired-DEG_analysis/", z, "_", gmt_db, "_GSEA.rds"),"\n")
  
  gmtLR = lapply(gmtList, pullDFrames, list_level=2)
  gmtLR = setNames(do.call(c, gmtLR), unlist(sapply(gmtLR, names)))
  
  gmtPair = lapply(gmtList, pullDFrames, list_level=3)
  gmtPair = setNames(do.call(c, gmtPair), unlist(sapply(gmtPair, names)))
  
  # heatmap of terms repeated across neuronal groups
  col_annot = cbind.data.frame("annotation"=factor(sapply(strsplit(names(gmtLR), "_"), function(x) x[[1]]), 
                                                   levels=c("sm","se"), labels=c("PRECAST","Seurat")),
                               "match_group"=c(rep("L2_L3_L4", 4), rep("L5", 2), rep("L6", 2)))
  rownames(col_annot) = names(gmtLR)
  annot_colors = list("annotation"=c("PRECAST"="dodgerblue4", "Seurat"="orange"),
                      "match_group"=c("L2_L3_L4"="#8ec7b7", "L5"=cpList$smoothed.light[["L5"]], 
                                      "L6"=cpList$smoothed.light[["L6"]]))
  
  LR.terms = do.call(rbind, lapply(names(gmtLR), function(x) {
    mutate(gmtLR[[x]], gene_set=x) %>% select(gene_set, pathway, NES)
  })
  ) %>% tidyr::pivot_wider(names_from="gene_set", values_from="NES", values_fill=0)
  LR.terms.m = as.matrix(LR.terms[,-1])
  #set and fix rownames
  rownames(LR.terms.m) = LR.terms$pathway
  if(gmt_db=="Reactome") {
    p.names = sapply(strsplit(LR.terms$pathway, " R-HSA-"), function(x) x[[1]])
  } 
  if(gmt_db=="WikiPathways") {
    p.names = sapply(strsplit(LR.terms$pathway, " WP"), function(x) x[[1]])
  }
  #p.names = sapply(p.names, function(x) paste(strwrap(x, width=75), collapse="\n"))
  p.names = sapply(p.names, function(x) {
    c1 = unlist(strwrap(x, width=75))
    if(length(c1)>1) {
      return(paste(c1[[1]],"[...]"))
    } else {
      return(c1[[1]])
    }
  })
  rownames(LR.terms.m) <- p.names
  #find and place missing columns
  missing.cols = setdiff(names(gmtLR), colnames(LR.terms.m))
  if(length(missing.cols)>0) {
    mm = matrix(0, nrow=nrow(LR.terms.m), ncol=length(missing.cols))
    colnames(mm) = missing.cols
    LR.terms.m = cbind(LR.terms.m, mm)
  }
  color.limits = ceiling(max(abs(LR.terms.m)))
  if(color.limits-max(abs(LR.terms.m))>.5) color.limits = color.limits-.5
  seq.len = color.limits*4
  pal.len = seq.len-1
  phm = pheatmap(LR.terms.m, color=colorRampPalette(rev(RColorBrewer::brewer.pal(n=5, "RdBu")))(pal.len),
                 breaks= seq(-color.limits, color.limits, length.out=seq.len),
                 fontsize_row = 7, treeheight_col = 12, treeheight_row = 12,
                 annotation_col= col_annot, annotation_colors= annot_colors,
                 clustering_method="complete", silent=T, angle_col=45,
                 main="NES of sig. terms as fill color (0 if not sig.)")
  
  # find JC of terms within L-R result based on the genes comprising the leading edge of the term
  # those results will be filtered to the pathways that also have a LR DEG
  jc.df = pullJC(gmtList, "LR")
  group_by(jc.df, orig_group) %>% tally() #used to check output
  
  # define group order to be consistent with L2_L3_L4, L5, L6
  og = c("sm_L2","sm_L3.4","se_L2.3","se_L4","sm_L5","se_L5","sm_L6","se_L6")
  # subset group order based on groups with results
  og = intersect(og, unique(jc.df$orig_group))
  igraphList <- lapply(og, generateIGRAPH, .gmt_db=gmt_db)
  names(igraphList) <- og
  
  pdf(file=paste0("plots/07-2_L-R_paired-DEG_analysis/LR_fgsea_", gmt_db, "_", gsub("_", "-", gsub("\\.", "-", z)), "_pathway-igraph.pdf"), 
      height=8, width=8)
  grid.arrange(phm[[4]], top=paste(gsub("_"," ", z), gmt_db))
  for(i in names(igraphList)) {
    plot(simplify(igraphList[[i]][["igraph"]]), layout=igraphList[[i]][["layout"]], 
         edge.width=E(igraphList[[i]][["igraph"]])$jc*5, 
         vertex.label.family="sans", #vertex.label.font=2, 
         vertex.frame.color=NA, main=paste(paste0(gsub("_", " ", z),":"), igraphList[[i]][["title_text"]]),
         sub=igraphList[[i]][["sub_text"]])
    legend(lpList[[z]][[i]], legend=c("Depleted (with LR DEG)","Depleted (with LR paired DEG)",
                                  "Enriched (with LR DEG)","Enriched (with LR paired DEG)"),
           pch=16, pt.cex=1, cex=.7,
           col=c("#CFEBF7","skyblue","#FFC0B5","tomato"))
  }
  dev.off()
  cat("\nSaved pathway igraph to:", paste0("plots/07-2_L-R_paired-DEG_analysis/LR_fgsea_", gmt_db, "_", gsub("_", "-", gsub("\\.", "-", z)), "_pathway-igraph.pdf"), "\n")
}

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
