setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(pheatmap)
	library(gridExtra)
	library(grid)
	library(gtable)
	library(gridtext)
})

set.seed(123)

cpList = readRDS("plots/colorPalettes.rds")

# PRECAST smoothed
restr.results_sm <- read.csv("processed-data/07_dx_DE/layer-restricted-age_smoothed-k9-1663_compiled-results.csv", row.names=1) %>%
  mutate(sex= factor(sex, levels=c("F","M")),
         group= factor(group, levels=c("NTC.MDD","NTC.BPD","MDD.BPD")),
         smoothed=factor(cluster, levels=c("L1","L2","L3.4","L5","L6","WM")),
         dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased", "increased")))

## filter to sig thresholds
restr_sm_tmp = bind_rows(filter(restr.results_sm, cluster!="WM", adj.P.Val<.05),
                         filter(restr.results_sm, cluster=="WM", sex=="M", adj.P.Val<.05),
                         filter(restr.results_sm, cluster=="WM", sex=="F", adj.P.Val<.01))
# Seurat label transfer
restr.results_se <- read.csv("processed-data/07_dx_DE/layer-restricted-age_seurat-pc30-no-lowUMI_compiled-results.csv", row.names=1) %>%
  mutate(sex= factor(sex, levels=c("F","M")),
         group= factor(group, levels=c("NTC.MDD","NTC.BPD","MDD.BPD")),
         seurat_label_f=factor(cluster, levels=c("Micro.Vasc","Astro","L2.3","L4","L5","L6","Oligo","Inhb"),
                               labels=c("M.V","Astro","L2.3","L4","L5","L6","Oligo","Inhb")),
         dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased", "increased")))

## filter to sig thresholds
restr_se_tmp = bind_rows(filter(restr.results_se, cluster!="Oligo", adj.P.Val<.05),
                         filter(restr.results_se, cluster=="Oligo", sex=="M", adj.P.Val<.05),
                         filter(restr.results_se, cluster=="Oligo", sex=="F", adj.P.Val<.01))

# correlation heatmap to find matching groups
gene_set = intersect(restr.results_sm$gene_id, restr.results_se$gene_id)
length(gene_set) #13003

coefList = c("NTC.MDD_F","NTC.MDD_M","NTC.BPD_F","NTC.BPD_M","MDD.BPD_F","MDD.BPD_M")
names(coefList) = coefList

### t stat
plist <- lapply(coefList, function(x) {
  tmp = unlist(strsplit(x, "_"))
  target_group = tmp[[1]]
  target_sex = tmp[[2]]
  
  t1 = filter(restr.results_sm, group==target_group, sex==target_sex, 
              gene_id %in% gene_set) %>% 
    select(t, gene_id, smoothed) %>% 
    tidyr::pivot_wider(names_from="smoothed", values_from="t")
  m1 = as.matrix(t1[,-1])
  rownames(m1) = t1$gene_id
  
  t2 = filter(restr.results_se, group==target_group, sex==target_sex,  gene_id %in% gene_set) %>% 
    select(t, gene_id, seurat_label_f) %>% 
    tidyr::pivot_wider(names_from="seurat_label_f", values_from="t")
  m2 = as.matrix(t2[,-1])
  rownames(m2) = t2$gene_id
  
  cor.m = cor(m1, m2[rownames(m1),])
  
  phm = pheatmap(cor.m[c("L1","L2","L3.4","L5","L6","WM"),
                       c("M.V","Inhb","Astro","L2.3","L4","L5","L6","Oligo")], 
                 cluster_rows=F, cluster_cols=F, angle_col=0, silent=T,
                 main=paste(target_group, target_sex), fontsize=7)
  return(phm[[4]])
})

### text for pdf
pointList = list("L5, L6, and WM/Oligo have extremely similar L-R results from the two annotations, consistent with the high degree of spot-level overlap between each pair of annotations.",
                 "There is a low degree of overlap between PRECAST L1 L-R results and Seurat label Micro.Vasc or Astro L-R results, despite some overlap in spot-level annotations.",
                 "There is no clean way to look at the supra-granular comparisons between the PRECAST and Seurat label annotations. Both from the biological justification of overlapping importance of each of L2 vs L2.3 and L3.4 vs L2.3 and L3.4 vs L4, as well as from looking at the correlation of the t-statstics of the L-R results.",
                 "Altogether, we will look at a combined group of L2,L3,L4 for L-R DEG consistency, as well as the 3 separate groups of L5, L6, and WM.")
tlist = lapply(pointList, function(x) {
  textbox_grob(
    text = x,
    x = unit(0.5, "npc"),
    y = unit(0.5, "npc"),
    width = unit(0.8, "npc"), # Text will wrap to fit 70% of the viewport width
    hjust = 0.5,
    vjust = 0.5,
    gp = gpar(fontsize = 10, lineheight = 1))
})


# generate DEG lists
nlist = as.character(sapply(coefList, function(x) paste0(x, c("_dn","_up"))))
names(nlist) = nlist

sm_match = c("L1"="L1", "L2_L3_L4"="L2", "L2_L3_L4"="L3.4", "L5"="L5", "L6"="L6", "WM"="WM")
se_match = c("M.V"="M.V", "Astro"="Astro", "L2_L3_L4"="L2.3", "L2_L3_L4"="L4", "L5"="L5", "L6"="L6", "Inhb"="Inhb", "WM"="Oligo")


pairList = lapply(nlist, function(x) {
  tmp = unlist(strsplit(x, "_"))
  target_group = tmp[[1]]
  target_sex = tmp[[2]]
  target_dir = tmp[[3]]
  target_dir2 = ifelse(target_dir=="dn", "decreased", "increased")
  
  outList = list()
  #LR_sig_any
  ## smoothed
  test1_sm = filter(restr_sm_tmp, group==target_group, sex==target_sex, dir==target_dir2)
  out1 = lapply(sm_match, function(x) {
    filter(test1_sm, smoothed== x)$gene_name
  })
  names(out1) = sm_match
  
  ## seurat
  test1_se = filter(restr_se_tmp, group==target_group, sex==target_sex, dir==target_dir2)
  out2 = lapply(se_match, function(x) {
    filter(test1_se, seurat_label_f== x)$gene_name
  })
  names(out2) = se_match
  
  ## concat
  outList$LR_sig_any = list(smoothed = out1, seurat= out2)
  
  #paired LR sig
  match_names = c("L2_L3_L4","L5","L6","WM")
  out3 = list()

  for(i in match_names) {
    orig_names = sm_match[grep(i, names(sm_match))]
    compare_names = se_match[grep(i, names(se_match))]
    if(length(orig_names)>1) {
      orig_genes =  unique(unlist(out1[orig_names]))
      compare_genes = unique(unlist(out2[compare_names]))
      out3[[i]] = intersect(orig_genes, compare_genes)
    } else {
      out3[[i]] = intersect(out1[[orig_names]], out2[[compare_names]])
    }
  }
  
  outList$LR_sig_paired = out3
  
  #solo LR sig
  ## smoothed
  out4 = list()
  for(i in match_names) {
    orig_names = sm_match[grep(i, names(sm_match))]
    if(length(orig_names)>1) {
      for(j in orig_names) {
        out4[[j]] = setdiff(out1[[j]], out3[[i]])
      }
    } else {
      out4[[orig_names]] = setdiff(out1[[orig_names]], out3[[i]])
    }
  }
  outList$LR_sig_solo$smoothed = out4
  
  ## seurat
  out5 = list()
  for(i in match_names) {
    orig_names = se_match[grep(i, names(se_match))]
    if(length(orig_names)>1) {
      for(j in orig_names) {
        out5[[j]] = setdiff(out2[[j]], out3[[i]])
      }
    } else {
      out5[[orig_names]] = setdiff(out2[[orig_names]], out3[[i]])
    }
  }
  outList$LR_sig_solo$seurat = out5

  return(outList)
})

saveRDS(pairList, file="processed-data/07-2_L-R_paired-DEG_analysis/LR-paired_lists.rds")
cat("\n\nSaved paired DEG list to: processed-data/07-2_L-R_paired-DEG_analysis/LR-paired_lists.rds\n")



# make bar plot of paired genes
bar.df = do.call(rbind, lapply(names(pairList), function(x) {
  tmp = unlist(strsplit(x, "_"))
  target_group = tmp[[1]]
  target_sex = tmp[[2]]
  target_dir = tmp[[3]]
  
  tmplist = names(pairList[[x]][["LR_sig_paired"]])
  tmp.df = do.call(rbind, lapply(tmplist, function(y) {
    p_length = length(pairList[[x]][["LR_sig_paired"]][[y]]) 
    y.sm = sm_match[grep(y, names(sm_match))]
    if(length(y.sm)>1) {
      sm_length = length(unique(unlist(pairList[[x]][["LR_sig_solo"]]$smoothed[y.sm])))
    } else sm_length = length(pairList[[x]][["LR_sig_solo"]]$smoothed[[y.sm]])
    
    y.se = se_match[grep(y, names(se_match))]
    if(length(y.se)>1) {
      se_length = length(unique(unlist(pairList[[x]][["LR_sig_solo"]]$seurat[y.se])))
    } else se_length = length(pairList[[x]][["LR_sig_solo"]]$seurat[[y.se]]) #5
    
    #format as dataframe
    as.data.frame(list(match_group = rep(y, 3),
                       set_type = c("paired","sm only","se only"),
                       set_size = c(p_length, sm_length, se_length)))
  })) %>% mutate(coef= paste(target_group, target_sex),
                 dir= target_dir)
  
  return(tmp.df)
})) %>% mutate(coef = factor(coef, levels=c("NTC.MDD F","NTC.MDD M","NTC.BPD F","NTC.BPD M","MDD.BPD F","MDD.BPD M")),
               match_group = factor(match_group, levels=c("L2_L3_L4","L5","L6","WM")),
               set_type = factor(set_type, levels=c("se only","paired","sm only"),
                                 labels=c("Seurat\nonly","Paired","PRECAST\nonly")))

p1 <- ggplot(bar.df, aes(y=factor(coef, levels=rev(levels(coef))), fill=set_type))+
  geom_bar(data=filter(bar.df, dir=="dn"), aes(x= -set_size), stat="identity", position="dodge",
           color="black", linewidth=.3, width=.7)+
  #geom_text(data=filter(df4.summ, dir=="dn"), aes(x= -mod_prop2, label= set_size, fill=NULL))+
  geom_bar(data=filter(bar.df, dir=="up"), aes(x= set_size), stat="identity", position="dodge", 
           color="black", linewidth=.3, width=.7)+
  #geom_text(data=filter(df4.summ, dir=="up"), aes(x= mod_prop2, label= set_size, fill=NULL))+
  scale_fill_manual(values=c("PRECAST\nonly"="dodgerblue4","Paired"="black", "Seurat\nonly"="orange"))+
  facet_wrap(vars(match_group), scales="free_x")+
  labs(x="# DEGs", fill="", y="")+
  geom_vline(aes(xintercept=0))+
  theme_minimal()+theme(axis.text.y=element_text(size=10), strip.text=element_text(size=12))


pointList = list("Each facet contains the data for one of the 4 groups chosen for paired L-R analysis. The y-axis contains the dx*sex comparison being examined and the x-axis indicates the number of DEGs past our signficance threshold, with the positive and negative values indicating if the DEG was increased or decreased, respectively.",
                 "The fill color indicates if the DEG was significant in only one annotation model or if the DEG was significant in the same direction in both annotation models. Fill color also matches the schematic.",
                 "These plots show that there is often considerable agreement between the DEGs of the two annotation strategies. However, there are also plenty of DEGs that were only significant in one model that we likely do not want to discard entirely.")

tlist2 = lapply(pointList, function(x) {
  textbox_grob(
    text = x,
    x = unit(0.5, "npc"),
    y = unit(0.5, "npc"),
    width = unit(0.8, "npc"), # Text will wrap to fit 70% of the viewport width
    hjust = 0.5,
    vjust = 0.5,
    gp = gpar(fontsize = 10, lineheight = 1))
})


# neuronal DEG vs glial DEG overlap
tmp = do.call(rbind, lapply(names(pairList), function(x) {
  tmp = unlist(strsplit(x, "_"))
  target_group = tmp[[1]]
  target_sex = tmp[[2]]
  target_dir = tmp[[3]]
  
  nrn.degs = unique(c(unlist(pairList[[x]][["LR_sig_any"]][["smoothed"]][c("L2","L3.4","L5","L6")]),
                      unlist(pairList[[x]][["LR_sig_any"]][["seurat"]][c("L2.3","L4","L5","L6")]))
                    )
  m1 = matrix(NA, nrow=4, ncol=2, dimnames= list(c("nrn_total","L1","Astro","M.V"), c("prop","total")))
  m1["nrn_total",] = c(1, length(nrn.degs))
  for(i in c("L1","Astro","M.V")) {
    if(i=="L1") {
      n_set = length(intersect(nrn.degs, unlist(pairList[[x]][["LR_sig_any"]][["smoothed"]][[i]])))
    } else n_set = length(intersect(nrn.degs, unlist(pairList[[x]][["LR_sig_any"]][["seurat"]][[i]])))
    m1[i,] = c(n_set/length(nrn.degs), n_set)
  }
  m1[,1] = round(m1[,1], 3)
  mutate(tibble::rownames_to_column(as.data.frame(m1), var="match_group"), 
         group=target_group, sex= target_sex, dir=target_dir)
})
) %>% mutate(match_group= factor(match_group, levels=c("nrn_total","L1","Astro","M.V"),
                                      labels=c("L2-6","L1","Astro","M.V")),
             group=factor(group, levels=c("NTC.MDD","NTC.BPD","MDD.BPD")),
             label_pos= ifelse(dir=="dn", -prop-.05, prop+.05),
             label_pos2= ifelse(dir=="dn" & label_pos>0, 0, label_pos),
             label_pos2= ifelse(dir=="up" & label_pos2<0, 0, label_pos2))

p2 = ggplot(tmp, aes(x=match_group, y=prop))+
  geom_bar(data=filter(tmp, dir=="dn"), aes(y=-prop), stat="identity", fill="skyblue", width=.6, color="black", linewidth=.3)+
  geom_text(data=filter(tmp, dir=="dn"), aes(y=label_pos2, label=total), vjust=1, size=3)+
  geom_bar(data=filter(tmp, dir=="up"), stat="identity", fill="tomato", width=.6, color="black", linewidth=.3)+
  geom_text(data=filter(tmp, dir=="up"), aes(y=label_pos2, label=total), vjust=0, size=3)+
  facet_grid(cols=vars(sex), rows=vars(group))+
  scale_y_continuous(expand = c(.1,.1))+
  geom_hline(aes(yintercept=0))+
  geom_vline(aes(xintercept=1.5), lty=2, linewidth=.3)+
  labs(y="prop. of L2-6 DEGs overlapping with un-paired glial DEGs")+
  theme_bw()+theme(axis.title.x= element_blank(), strip.text=element_text(size=10))

pointList = list("This plot examines the overlap between any neuronal (L2-L6) DEGs from either annotation and the glial DEGs of L1 (PRECAST), Micro.Vasc (Seurat), and Astro (Seurat).",
                 "For each dx*sex comparison (facets), the y-axis indicates the proportion of L2-6 DEGs that were also DEGs in the x-axis groups. Positive/ negative values indicate whether DEGs were increased/ decreased, respectively. To help orient, the vertical dashed line separates the gene set being used in the denominator (on the left; L2-6 DEGs), versus the glial DEGs being used in the numerator (on the right).",
                 "The total number of DEGs in each bar is given by the text value, to help discern when a large proportion of L2-6 DEGs being identified as glial DEGs is due to low overall number of L2-6 DEGs. These values are also direction (increased/ decreased) specific."
                 )

tlist3 = lapply(pointList, function(x) {
  textbox_grob(
    text = x,
    x = unit(0.5, "npc"),
    y = unit(0.5, "npc"),
    width = unit(0.8, "npc"), # Text will wrap to fit 70% of the viewport width
    hjust = 0.5,
    vjust = 0.5,
    gp = gpar(fontsize = 10, lineheight = 1))
})


#DEG tables
tt <- ttheme_default(
  core = list(fg_params=list(cex = .6)),
  colhead = list(fg_params=list(cex = .6)),
  rowhead = list(fg_params=list(cex = .6))
  )

glist = lapply(names(pairList), function(x){
  tmp = unlist(strsplit(x, "_"))
  target_group = tmp[[1]]
  target_sex = tmp[[2]]
  target_dir = tmp[[3]]
  
  if(target_dir=="dn") {
    fill_color="#CFEBF7"
  } else fill_color="#FFC0B5"
  if(target_dir=="dn") {
    target_dir2="Decreased"
  } else target_dir2="Increased"
  
  
  #precast table
  fills = rep(fill_color, length(cpList$smoothed.light)*2)
  fills[2] = "white"
  tt1 = modifyList(tt, list(colhead= list(bg_params=list(fill=cpList$smoothed.bright)),
                            core= list(bg_params=list(fill=fills)
                            ))
  )
  
  df1 = rbind.data.frame(as.data.frame(lapply(pairList[[x]][["LR_sig_any"]]$smoothed, length)),
                         c(rep("", length(pairList[[x]][["LR_sig_any"]]$smoothed)))
  )
  rownames(df1) = c("PRECAST DEG", "PRECAST only")
  
  for(i in names(pairList[[x]][["LR_sig_solo"]]$smoothed)) {
    df1[2,i] = length(pairList[[x]][["LR_sig_solo"]]$smoothed[[i]])
  }
  
  gt1 = tableGrob(df1, theme = tt1)
  
  
  #seurat table
  color.mtx = rbind(rep(fill_color, length(cpList$transfer.light)),
                    rep(fill_color, length(cpList$transfer.light)),
                    cpList$transfer.bright[c("Inhb","Micro.Vasc","Astro","L2.3","L4","L5","L6","Oligo")])
  color.mtx[1,1:3] ="white"
  fills2 = as.character(color.mtx)
  
  tt2 = modifyList(tt, list(colhead= list(bg_params=list(fill=cpList$transfer.bright[c("Inhb","Micro.Vasc","Astro","L2.3","L4","L5","L6","Oligo")])),
                            core= list(bg_params=list(fill=fills2),
                                       fg_params=list(fontface=c("plain","plain","bold")))))
  
  df2 = rbind.data.frame(as.data.frame(lapply(pairList[[x]][["LR_sig_any"]]$seurat, length)),
                         c(rep("", length(pairList[[x]][["LR_sig_any"]]$seurat)))
  )
  
  for(i in names(pairList[[x]][["LR_sig_solo"]]$seurat)) {
    df2[2,i] = length(pairList[[x]][["LR_sig_solo"]]$seurat[[i]])
  }
  df2 = df2[c(2,1),c("Inhb","M.V","Astro","L2.3","L4","L5","L6","Oligo")]
  rownames(df2) = c("Seurat only","Seurat DEG")
  df2 = rbind.data.frame(df2, colnames(df2))
  rownames(df2)[3] =""
  
  gt2 = tableGrob(df2, cols=NULL, theme=tt2)
  
  
  #paired row
  tt3 = modifyList(tt, list(rowhead= list(fg_params= list(fontface=c("bold.italic"))),
                            core= list(bg_params= list(fill=c("#8ec7b7",cpList$smoothed.light[["L5"]], 
                                                              cpList$smoothed.light[["L6"]], 
                                                              cpList$smoothed.light[["WM"]])),
                                       fg_params= list(fontface=c("bold.italic"))))
  )
  df3 = t(as.data.frame(sapply(pairList[[x]][["LR_sig_paired"]], length)))
  gt3 = tableGrob(df3, rows="Paired DEG", cols=NULL, theme= tt3)
  
  #combine with seurat
  gt4 = gtable_combine(gt3, gt2, along=2)
  
  gt4$layout[c(3,7),"l"] = 5
  gt4$layout[c(3,7),"r"] = 6
  gt4$layout[c(4:6,8:10),"l"] = gt4$layout[c(4:6,8:10),"l"]+4
  gt4$layout[c(4:6,8:10),"r"] = gt4$layout[c(4:6,8:10),"r"]+4
  
  
  #combine with precast
  gt5 = gtable_combine(gt1, gt4, along=2)
  
  gt5$layout[7:42,"l"] = gt5$layout[7:42,"l"]+2
  gt5$layout[7:42,"r"] = gt5$layout[7:42,"r"]+2
  
  #return(gt5)
  #ggsave(file="plots/07-2_L-R_paired-DEG_analysis/test.png", gt5)
  
  #add direction as title
  title <- textGrob(target_dir2, gp=gpar(fontsize=10))
  padding <- unit(5,"mm")
  table <- gtable_add_rows(
    gt5, 
    heights = grobHeight(title) + padding,
    pos = 0)
  
  table <- gtable_add_grob(
    table, 
    title, 
    1, 1, 1, ncol(table))
  table$layout
  table$layout[107,"r"] = 2
  
  return(table)
})


#save summary PDF
pdf(file="plots/07-2_L-R_paired-DEG_analysis/LR-paired_summary.pdf", width=8, height=6)
grid.arrange(arrangeGrob(grobs=plist, ncol=2), 
             arrangeGrob(grobs=tlist, layout_matrix=matrix(c(1,2,3,3,4))), 
             layout_matrix=matrix(c(1,1,2), ncol=3), top="t statistic correlation between L-R annotations")
grid.arrange(p1,
             arrangeGrob(grobs=tlist2, layout_matrix=matrix(c(1,1,1,2,2,3,3))),
             layout_matrix=matrix(c(1,1,2), ncol=3), top="Number of paired L-R DEGs across matching groups")
grid.arrange(p2,
             arrangeGrob(grobs=tlist3, layout_matrix=matrix(c(1,1,2,2,2,2,2,3,3,3))),
             layout_matrix=matrix(c(1,1,2), ncol=3))
grid.arrange(arrangeGrob(grobs=glist[1:2], ncol=2, top="NTC.MDD F"),
             arrangeGrob(grobs=glist[3:4], ncol=2, top="NTC.MDD M"), ncol=1)
grid.arrange(arrangeGrob(grobs=glist[5:6], ncol=2, top="NTC.BPD F"),
             arrangeGrob(grobs=glist[7:8], ncol=2, top="NTC.BPD M"), ncol=1)
grid.arrange(arrangeGrob(grobs=glist[9:10], ncol=2, top="MDD.BPD F"),
             arrangeGrob(grobs=glist[11:12], ncol=2, top="MDD.BPD M"), ncol=1)
dev.off()
cat("\n\nSaved summary PDF to: plots/07-2_L-R_paired-DEG_analysis/LR-paired_summary.pdf\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
