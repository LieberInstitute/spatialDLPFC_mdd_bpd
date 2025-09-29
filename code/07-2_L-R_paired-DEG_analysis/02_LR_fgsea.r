setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(pheatmap)
	library(fgsea)
	library(enrichR)
	library(gridExtra)
	library(ggrastr)
	library(grid)
	library(gtable)
	library(gridtext)
})

set.seed(123)

setEnrichrSite("Enrichr") # Human genes
dbs <- listEnrichrDbs()

source("code/07-2_L-R_paired-DEG_analysis/fgsea_functions.r")
source("code/07-2_L-R_paired-DEG_analysis/manhattan_functions.r")
cpList = readRDS("plots/colorPalettes.rds")

# load L-R DEG lists
pairList = readRDS("processed-data/07-2_L-R_paired-DEG_analysis/LR-paired_logFC-0.3_lists.rds")

sm_match = c("L1"="L1", "L2_L3_L4"="L2", "L2_L3_L4"="L3.4", "L5"="L5", "L6"="L6", "WM"="WM")
se_match = c("M.V"="M.V", "Astro"="Astro", "L2_L3_L4"="L2.3", "L2_L3_L4"="L4", "L5"="L5", "L6"="L6", "Inhb"="Inhb", "WM"="Oligo")

# load L-R results for logFC ranking
restr.results_sm <- read.csv("processed-data/07_dx_DE/layer-restricted-age_smoothed-k9-1663_compiled-results.csv", row.names=1) %>%
  mutate(sex= factor(sex, levels=c("F","M")),
         group= factor(group, levels=c("NTC.MDD","NTC.BPD","MDD.BPD")),
         smoothed=factor(cluster, levels=c("L1","L2","L3.4","L5","L6","WM")),
         dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased", "increased")))

restr.results_se <- read.csv("processed-data/07_dx_DE/layer-restricted-age_seurat-pc30-no-lowUMI_compiled-results.csv", row.names=1) %>%
  mutate(sex= factor(sex, levels=c("F","M")),
         group= factor(group, levels=c("NTC.MDD","NTC.BPD","MDD.BPD")),
         seurat_label_f=factor(cluster, levels=c("Micro.Vasc","Astro","L2.3","L4","L5","L6","Oligo","Inhb"),
                               labels=c("M.V","Astro","L2.3","L4","L5","L6","Oligo","Inhb")),
         dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased", "increased")))


# what if looked at la lr overlap genes too
saveList <- readRDS("processed-data/07_dx_DE/LA-LR-overlap_lists.rds")


# loop for all dx*sex comparisons and match groups
coefList = c("NTC.MDD_F","NTC.MDD_M","NTC.BPD_F","NTC.BPD_M","MDD.BPD_F","MDD.BPD_M")
names(coefList) = coefList

matchList = c("L2_L3_L4","L5","L6")
names(matchList) = matchList

#gmt_db = "Reactome"
gmt_db = "WikiPathways"

if(gmt_db=="Reactome") gmt = .read_gmt("Reactome_2022")
if(gmt_db=="WikiPathways") gmt = .read_gmt("WikiPathways_2024_Human")

gmt.list = group_by(gmt, term) %>% summarise(gene=list(gene)) %>%
	tibble::deframe()


## prepare outlist for gtables
glist <- list()
glist2 <- list() #for LA.LR both annotations tables
tt <- ttheme_default(
  core = list(fg_params=list(cex = .6)),
  colhead = list(fg_params=list(cex = .6)),
  rowhead = list(fg_params=list(cex = .6))
)


## execute fgsea and save, generate tables of # of results
for(z in coefList) {
  cat("\n************\n")
  cat(z,"\n")

  tmp = unlist(strsplit(z, "_"))
  target_group = tmp[[1]]
  target_sex = tmp[[2]]
  
  # skip generating if already exists
  if(!paste0(z,"_", gmt_db, "_GSEA.rds") %in% list.files("processed-data/07-2_L-R_paired-DEG_analysis/")) {
	gmtList = lapply(matchList, generateFGSEA, .gmtList=gmt.list)
	saveRDS(gmtList, paste0("processed-data/07-2_L-R_paired-DEG_analysis/", z, "_", gmt_db, "_GSEA.rds"))
	cat("\nSaved", gmt_db, "results to:", paste0("processed-data/07-2_L-R_paired-DEG_analysis/", z, "_", gmt_db, "_GSEA.rds"),"\n")
  } else {
	gmtList <- readRDS(paste0("processed-data/07-2_L-R_paired-DEG_analysis/", z, "_", gmt_db, "_GSEA.rds"))
	cat("\nLoad saved", gmt_db, "results from:", paste0("processed-data/07-2_L-R_paired-DEG_analysis/", z, "_", gmt_db, "_GSEA.rds"),"\n")
  }
  gmtPadj = lapply(gmtList, pullDFrames, list_level=1)
  gmtPadj = setNames(do.call(c, gmtPadj), unlist(sapply(gmtPadj, names)))
  
  gmtLR = lapply(gmtList, pullDFrames, list_level=2)
  gmtLR = setNames(do.call(c, gmtLR), unlist(sapply(gmtLR, names)))
  
  gmtPair = lapply(gmtList, pullDFrames, list_level=3)
  gmtPair = setNames(do.call(c, gmtPair), unlist(sapply(gmtPair, names)))

#  # skip generating if already exists
#  if(!paste0(z,"_WikiPathways_GSEA.rds") %in% list.files("processed-data/07-2_L-R_paired-DEG_analysis/")) {
#	wikiList = lapply(matchList, generateFGSEA, .gmtList=wiki.list)
#	saveRDS(wikiList, paste0("processed-data/07-2_L-R_paired-DEG_analysis/", z, "_WikiPathways_GSEA.rds"))
#	cat("\nSaved WikiPathways results to:", paste0("processed-data/07-2_L-R_paired-DEG_analysis/", z, "_WikiPathways_GSEA.rds"),"\n")
#  } else {
#	wikiList <- readRDS(paste0("processed-data/07-2_L-R_paired-DEG_analysis/", z, "_WikiPathways_GSEA.rds"))
#	cat("\nLoad saved wikipathways results from:", paste0("processed-data/07-2_L-R_paired-DEG_analysis/", z, "_WikiPathways_GSEA.rds"),"\n")
#  }
#  wikiPadj = lapply(wikiList, pullDFrames, list_level=1)
#  wikiPadj = setNames(do.call(c, wikiPadj), unlist(sapply(wikiPadj, names)))
#
#  wikiLR = lapply(wikiList, pullDFrames, list_level=2)
#  wikiLR = setNames(do.call(c, wikiLR), unlist(sapply(wikiLR, names)))
#
#  wikiPair = lapply(wikiList, pullDFrames, list_level=3)
#  wikiPair = setNames(do.call(c, wikiPair), unlist(sapply(wikiPair, names)))

  ## gtables ----
  start.df = cbind.data.frame("n_sig_terms_padj.05"=sapply(gmtPadj, nrow),
                              "n_sig_terms_with_LR.DEG"=sapply(gmtLR, nrow),
                              "n_sig_terms_with_paired.LR.DEG"=sapply(names(gmtLR), function(x) {
                                length(intersect(gmtLR[[x]]$pathway,gmtPair[[x]]$pathway))
                              }))
  
#  start.df = cbind.data.frame("n_sig_terms_padj.05"=sapply(wikiPadj, nrow),
#                              "n_sig_terms_with_LR.DEG"=sapply(wikiLR, nrow),
#                              "n_sig_terms_with_paired.LR.DEG"=sapply(names(wikiLR), function(x) {
#                                length(intersect(wikiLR[[x]]$pathway,wikiPair[[x]]$pathway))
#                              }))

  rnames = rownames(start.df)
  deg.df = do.call(rbind, lapply(rnames, function(i) {
    tmp = unlist(strsplit(i, "_"))
    map_annot = as.character(factor(tmp[[1]], levels=c("sm","se"), labels=c("smoothed","seurat")))
    map_orig = tmp[[2]]
    map_match = if(map_annot=="smoothed") {
      names(sm_match)[sm_match==map_orig]
    } else {
      names(se_match)[se_match==map_orig]
    }
    tibble::enframe(list(n.LR.degs = length(union(pairList[[paste0(z,"_dn")]][[1]][[map_annot]][[map_orig]],
                                                  pairList[[paste0(z,"_up")]][[1]][[map_annot]][[map_orig]])),
                         n.paired.LR.degs = length(union(pairList[[paste0(z,"_dn")]][[2]][[map_match]],
                                                         pairList[[paste0(z,"_up")]][[2]][[map_match]]))
    )) %>%
      mutate(value=as.numeric(value), orig_group=i)
  })) %>% tidyr::pivot_wider(names_from="name", values_from="value") %>%
    tibble::column_to_rownames(var="orig_group")
  
  t1 = cbind(start.df, deg.df)[,c(1,4,2,5,3)]
  t1.fill = matrix(c(rep("grey95", ncol(t1)), rep("grey90", ncol(t1))), nrow=nrow(t1), ncol=ncol(t1), byrow=T)

  t2 = do.call(rbind, lapply(rownames(t1), function(i) {
    tmp = unlist(strsplit(i, "_"))
    map_annot = as.character(factor(tmp[[1]], levels=c("sm","se"), labels=c("smoothed","seurat")))
    map_orig = tmp[[2]]
    map_match = if(map_annot=="smoothed") {
      names(sm_match)[sm_match==map_orig]
    } else {
      names(se_match)[se_match==map_orig]
    }
    data.frame("annot"=map_annot, "match"=map_match, "orig_group"=map_orig)
  })
  )
  t2.fill = cbind(c("#87A6C4","#87A6C4","#FFD27F","#FFD27F","#87A6C4","#FFD27F","#87A6C4","#FFD27F"),
                  c(rep("#8ec7b7", 4), rep(cpList$smoothed.light[["L5"]], 2), rep(cpList$smoothed.light[["L6"]], 2)),
                  c(cpList$smoothed.light[c("L2","L3.4")], cpList$transfer.light[c("L2.3","L4")], 
                    cpList$smoothed.light[["L5"]], cpList$transfer.light[["L5"]],
                    cpList$smoothed.light[["L6"]], cpList$transfer.light[["L6"]]))

  fills = as.character(cbind(t2.fill, t1.fill))
  tt1 = modifyList(tt, list(core= list(bg_params=list(fill=fills))))

  g = tableGrob(cbind(t2, t1), cols=c("annotation","neuron group","source","# sig. terms",
                                      "# LR DEGs","# sig. terms\nwith LR DEG",
                                      "# paired LR degs","# sig. terms\nwith paired LR DEG"),
                rows=NULL, theme=tt1)
  glist[[z]] = gtable_add_grob(g, 
                               grobs = rectGrob(gp=gpar(fill=NA, 
                                                        lwd=2)), 
                               t = 2, b = nrow(g), l = 6, r = 6)

  ### add second table with LA.LR both annotations genes
  la.lr.genes = c(saveList$LA.LR_both.annotations[[paste(target_group, target_sex, "dn", sep="_")]],
                  saveList$LA.LR_both.annotations[[paste(target_group, target_sex, "up", sep="_")]])
  LA.LR.terms = lapply(gmtPadj, function(x) {
    tidyr::separate_rows(x, leadingEdge2, sep="/") %>% 
      filter(leadingEdge2 %in% la.lr.genes) %>%
      pull(pathway) %>% unique()
  })
  start.df2 = cbind.data.frame("n_sig_terms_padj.05"=sapply(gmtPadj, nrow),
                               "n_LA.LR.both_DEGs"=length(la.lr.genes),
                               "n_sig_terms_with_LA.LR.both"=sapply(LA.LR.terms, length)
                              )
  fills2 = as.character(cbind(t2.fill, t1.fill[,1:3]))
  tt2 = modifyList(tt, list(core= list(bg_params=list(fill=fills2))))
  glist2[[z]] = tableGrob(cbind(t2, start.df2), cols=c("annotation","neuron group","source","# sig. terms",
                                  "# LA & LR DEGs","# sig. terms with\n LA & LR DEG"),
            rows=NULL, theme=tt2)
}

# example manhattan plots
# reactome
if(gmt_db=="Reactome") {
	cat("\n\n\nExample manhattan plots using NTC.MDD F Seurat L4\n") #for reactome
	## read in results
	gmtList <- readRDS("processed-data/07-2_L-R_paired-DEG_analysis/NTC.MDD_F_Reactome_GSEA.rds")

	gmtPadj = lapply(gmtList, pullDFrames, list_level=1)
	gmtPadj = setNames(do.call(c, gmtPadj), unlist(sapply(gmtPadj, names)))

	padj.only = "Neurotransmitter Receptors And Postsynaptic Signal Transmission R-HSA-112314"
	lr.only = "Neuronal System R-HSA-112316"
	pair.ex = "Cytokine Signaling In Immune System R-HSA-1280215"

	gsea.plot = formatManhattan("NTC.MDD", "F", "L4", "seurat", gmt.list, gmtPadj,
                            gsea_terms = c(padj.only, lr.only, pair.ex))

	v.df = filter(restr.results_se, group=="NTC.MDD", sex=="F", seurat_label_f=="L4")

	tx3.text = "Volcano plot (right): In most cases, filtering by paired L-R DEGs is more restrictive (generates fewer results). However, for the individual PRECAST clusters/ Seurat labels that comprise L2_L3_L4 like this example, our approach of looking at the union of L-R DEGs within a single annotation means that some paired L-R DEGs (black asterisk) are not statistically significant in the group for which the fgsea was performed (see black asterisks below the dashed red line)."
}
# wikipathways
if(gmt_db=="WikiPathways") {
	cat("\n\n\nExample manhattan plots using NTC.MDD F Seurat L5\n")
	# read in results
	gmtList <- readRDS("processed-data/07-2_L-R_paired-DEG_analysis/NTC.MDD_F_WikiPathways_GSEA.rds")

	gmtPadj = lapply(gmtList, pullDFrames, list_level=1)
	gmtPadj = setNames(do.call(c, gmtPadj), unlist(sapply(gmtPadj, names)))

	padj.only = "Complement System In Neuronal Development And Plasticity WP5090"
	lr.only = "Microglia Pathogen Phagocytosis Pathway WP3937"
	pair.ex = "Macrophage Stimulating Protein MSP Signaling WP5353"

	gsea.plot = formatManhattan("NTC.MDD", "F", "L5", "seurat", gmt.list, gmtPadj,
                            gsea_terms = c(padj.only, lr.only, pair.ex))

	v.df = filter(restr.results_se, group=="NTC.MDD", sex=="F", seurat_label_f=="L5")

	tx3.text = "Volcano plot (right): As this example shows, it is possible for a L-R DEG to be significantly different in the opposite direction of the GSEA result."
}


m1 = ggplot(filter(gsea.plot[["plotting"]], present==T, term==padj.only))+
  geom_segment(aes(x=rank, xend=rank, y=0, yend=weight), linewidth=.3, color="grey")+
  geom_text(data=filter(gsea.plot[["plotting"]], present2!="NS", term==padj.only),
            aes(x=rank, y=weight, label="*", color=present2), size=6)+
  scale_x_continuous(limits=c(0,max(gsea.plot[["plotting"]]$rank)), expand=c(.01,.01))+
  scale_color_manual(values=c("LR"="red3", "paired"="black"), guide="none")+
  geom_segment(data=data.frame(x1=0, x2=max(gsea.plot[["plotting"]]$rank), y1=0, y2=0),
               aes(x=x1, xend=x2, y=y1, yend=y2), color="grey")+
  geom_text(data=filter(gsea.plot[["labels"]], term==padj.only), aes(x=median(gsea.plot[["plotting"]]$rank), y=.1, label=y_label), 
            vjust=0, size=3)+
  theme_minimal()+labs(title=padj.only, y="logFC")+
  theme(axis.title=element_text(size=8), axis.text=element_text(size=7),
        plot.title=element_text(size=10, face="bold"))

v1 = inner_join(v.df, select(filter(gsea.plot[["plotting"]], term==padj.only), id, weight, present),
               by=c("gene_name"="id","logFC"="weight")) %>%
  mutate(point_col= factor(present, levels=c(FALSE,TRUE), labels=c("not present", "in gene set\nfor term")))

p1 <- ggplot(v1, aes(x=logFC, y=-log10(adj.P.Val), color=point_col))+
  geom_point(data=filter(v1, present==F), size=.5)+
  geom_point(data=filter(v1, present==T), size=.5)+
  geom_hline(aes(yintercept=-log10(.05)), lty=2, color="red3")+
  geom_text(data=data.frame("x1"=2.7, "y1"=-log10(.05), "lab1"="adj. p<.05"),
            aes(x=x1, y=y1, label=lab1), color="red3", size=3, fontface="italic", 
            vjust=0, nudge_y = .05)+
  scale_color_manual("Gene status:", values=c("grey80","grey50"),
                     guide=guide_legend(override.aes = list(size=2)))+
  labs(title="Volcano showing genes contributing to GSEA result")+
  theme_bw()+theme(text=element_text(size=8), legend.box.spacing= unit(1,"pt"))

m2 = ggplot(filter(gsea.plot[["plotting"]], present==T, term==lr.only))+
  geom_segment(aes(x=rank, xend=rank, y=0, yend=weight), linewidth=.3, color="grey")+
  geom_text(data=filter(gsea.plot[["plotting"]], present2!="NS", term==lr.only),
            aes(x=rank, y=weight, label="*", color=present2), size=6, vjust=1)+
  scale_x_continuous(limits=c(0,max(gsea.plot[["plotting"]]$rank)), expand=c(.01,.01))+
  scale_color_manual(values=c("LR"="red3", "paired"="black"), guide="none")+
  geom_segment(data=data.frame(x1=0, x2=max(gsea.plot[["plotting"]]$rank), y1=0, y2=0),
               aes(x=x1, xend=x2, y=y1, yend=y2), color="grey")+
  geom_text(data=filter(gsea.plot[["labels"]], term==lr.only), aes(x=median(gsea.plot[["plotting"]]$rank), y=.1, label=y_label), 
            vjust=0, size=3)+
  theme_minimal()+labs(title=lr.only, y="logFC")+
  theme(axis.title=element_text(size=8), axis.text=element_text(size=7),
        plot.title=element_text(size=10, face="bold"))

v2 = inner_join(v.df, select(filter(gsea.plot[["plotting"]], term==lr.only), id, weight, present, present2),
               by=c("gene_name"="id","logFC"="weight")) %>%
  mutate(point_col= factor(paste(present, present2), levels=c("FALSE NS", "TRUE NS", "TRUE LR"),
                           labels=c("not present","in gene set\nfor term","in gene set\n& L-R DEG")))

p2 <- ggplot(v2, aes(x=logFC, y=-log10(adj.P.Val), color=point_col))+
  geom_point(data=filter(v2, present==F), size=.5)+
  geom_point(data=filter(v2, present==T, present2=="NS"), size=.5)+
  geom_hline(aes(yintercept=-log10(.05)), lty=2, color="red3")+
  geom_text(data=data.frame("x1"=2.7, "y1"=-log10(.05), "lab1"="adj. p<.05"),
            aes(x=x1, y=y1, label=lab1), color="red3", size=3, fontface="italic", 
            vjust=0, nudge_y = .05)+
  geom_text(data=filter(v2, present==T, present2!="NS"), aes(label="*"), size=5, fontface="bold")+
  scale_color_manual("Gene status:", values=c("grey80","grey50","red3"),
                     guide=guide_legend(override.aes = list(label="*", size=c(2,2,5))))+
  labs(title="Volcano showing genes contributing to GSEA result")+
  theme_bw()+theme(text=element_text(size=8), legend.box.spacing= unit(1,"pt"),
                   legend.text= element_text(margin= margin(0,0,0,3,"pt")))

m3 = ggplot(filter(gsea.plot[["plotting"]], present==T, term==pair.ex))+
  geom_segment(aes(x=rank, xend=rank, y=0, yend=weight), linewidth=.3, color="grey")+
  geom_text(data=filter(gsea.plot[["plotting"]], present2!="NS", term==pair.ex),
            aes(x=rank, y=weight, label="*", color=present2), size=6)+
  scale_x_continuous(limits=c(0,max(gsea.plot[["plotting"]]$rank)), expand=c(.01,.01))+
  scale_y_continuous(expand=c(0,.25))+
  scale_color_manual(values=c("LR"="red3", "paired"="black"), guide="none")+
  geom_segment(data=data.frame(x1=0, x2=max(gsea.plot[["plotting"]]$rank), y1=0, y2=0),
               aes(x=x1, xend=x2, y=y1, yend=y2), color="grey")+
  geom_text(data=filter(gsea.plot[["labels"]], term==pair.ex), aes(x=median(gsea.plot[["plotting"]]$rank), y=.1, label=y_label), 
            vjust=0, size=3)+
  theme_minimal()+labs(title=pair.ex, y="logFC")+
  theme(axis.title=element_text(size=8), axis.text=element_text(size=7),
        plot.title=element_text(size=10, face="bold"))

v3 = inner_join(v.df, select(filter(gsea.plot[["plotting"]], term==pair.ex), id, weight, present, present2),
                by=c("gene_name"="id","logFC"="weight")) %>%
  mutate(point_col= factor(paste(present, present2), levels=c("FALSE NS", "TRUE NS", "TRUE LR", "TRUE paired"),
                           labels=c("not present","in gene set\nfor term","in gene set\n& L-R DEG","in gene set &\npaired L-R DEG")))

p3 <- ggplot(v3, aes(x=logFC, y=-log10(adj.P.Val), color=point_col))+
  geom_point(data=filter(v3, present==F), size=.5)+
  geom_point(data=filter(v3, present==T, present2=="NS"), size=.5)+
  geom_hline(aes(yintercept=-log10(.05)), lty=2, color="red3")+
  geom_text(data=data.frame("x1"=2.7, "y1"=-log10(.05), "lab1"="adj. p<.05"),
            aes(x=x1, y=y1, label=lab1), color="red3", size=3, fontface="italic", 
            vjust=0, nudge_y = .05)+
  geom_text(data=filter(v3, present==T, present2!="NS"), aes(label="*"), size=5, fontface="bold")+
  scale_color_manual("Gene status:", values=c("grey80","grey50","red3","black"),
                     guide=guide_legend(override.aes = list(label="*", size=c(2,2,5,5))))+
  labs(title="Volcano showing genes contributing to GSEA result")+
  theme_bw()+theme(text=element_text(size=8), legend.box.spacing= unit(1,"pt"),
                   legend.text= element_text(margin= margin(0,0,0,0,"pt")))


## text for manhattan examples
pointC = c("Since fgsea results are calculated based on the ranked logFC values, many fgsea results with padj<.05 do not contain any L-R DEGs. For our analysis, we have filtered out such results.<br>",
           "fgsea manhattan plot (top): Lowly ranked values (x-axis) have the highest, most positive logFC (y-axis).",
           "Volcano plot (right): Opposite from above, the logFC (x-axis) shows negative values on left and positive values on right.")

tx1 = textbox_grob(
  text = paste(pointC, collapse="<br><br>"),
  width = unit(0.8, "npc"), # Text will wrap to fit 80% of the viewport width
  hjust = 0.5, vjust = 0.5,
  halign = 0, valign=0,
  gp = gpar(fontsize = 10, lineheight = 1))


pointC2 = c("Significant (padj<.05) fgsea results were filtered to those that contained a L-R DEG in the leading edge genes. These L-R DEGs did not have to be paired (sig. in both annotations), allowing for unique PRECAST cluster or Seurat cell type label results to emerge (if present).<br>",
           "fgsea manhattan plot (top): The rank position of L-R DEGs is indicated with a red asterisk.",
           "Volcano plot (right): The L-R DEGs are indicated with a red asterisk.")

tx2 = textbox_grob(
  text = paste(pointC2, collapse="<br><br>"),
  width = unit(0.8, "npc"), # Text will wrap to fit 80% of the viewport width
  hjust = 0.5, vjust = 0.5,
  halign = 0, valign=0,
  gp = gpar(fontsize = 10, lineheight = 1))


pointC3 = c("Significant (padj<.05) fgsea results containing a L-R DEG were further stratified by those results that also contain a paired L-R DEG in the leading edge genes. This will help distinguish the smaller number of sig. results that only contained layer-specific L-R DEGs that were not reproducible in the other annotation strategy.<br>",
            "fgsea manhattan plot (top): The rank position of paired L-R DEGs is indicated with a black asterisk. Any additional layer-specific L-R DEGs are also indicated with a red asterisk.",
	tx3.text)
#            "Volcano plot (right): In most cases, filtering by paired L-R DEGs is more restrictive (generates fewer results). However, for the individual PRECAST clusters/ Seurat labels that comprise L2_L3_L4 like this example, our approach of looking at the union of L-R DEGs within a single annotation means that some paired L-R DEGs (black asterisk) are not statistically significant in the group for which the fgsea was performed (see black asterisks below the dashed red line).")
#	     "Volcano plot (right): As this example shows, it is possible for a L-R DEG to be significantly different in the opposite direction of the GSEA result.")

tx3 = textbox_grob(
  text = paste(pointC3, collapse="<br><br>"),
  width = unit(0.8, "npc"), # Text will wrap to fit 80% of the viewport width
  hjust = 0.5, vjust = 0.5,
  halign = 0, valign=0,
  gp = gpar(fontsize = 10, lineheight = 1))


# save summary pdf
pdf(file=paste0("plots/07-2_L-R_paired-DEG_analysis/LR_fgsea_", gmt_db, "_summary.pdf"), width=8, height=6)
grid.arrange(m1, rasterize(p1, dpi=150), tx1,
             layout_matrix=rbind(c(1,1,1,1),c(3,3,2,2),c(3,3,2,2)))
grid.arrange(m2, rasterize(p2, dpi=150), tx2, 
             layout_matrix=rbind(c(1,1,1,1),c(3,3,2,2),c(3,3,2,2)))
grid.arrange(m3, rasterize(p3, dpi=150), tx3,
             layout_matrix=rbind(c(1,1,1,1),c(3,3,2,2),c(3,3,2,2)))
for(i in names(glist)) {
  grid.arrange(gtable_combine(glist[[i]], glist2[[i]], along=2), top=paste(gsub("_", " ", i), gmt_db))
  #grid.arrange(glist[[i]], top=gsub("_", " ", i))
}
dev.off()
cat("\nSummary PDF saved to:", paste0("plots/07-2_L-R_paired-DEG_analysis/LR_fgsea_", gmt_db, "_summary.pdf"),"\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
