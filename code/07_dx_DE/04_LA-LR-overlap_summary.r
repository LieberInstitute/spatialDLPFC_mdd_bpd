setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	#library(SpatialExperiment)
	#library(edgeR)
	library(dplyr)
	library(ggplot2)
	#library(scater)
	#library(pheatmap)
	library(gridExtra)
	library(grid)
	library(gtable)
	library(ggrastr)
})

set.seed(123)

cpList = readRDS("plots/colorPalettes.rds")


# load DE model results
## PRECAST smoothed
adj.results_sm = read.csv("processed-data/07_dx_DE/layer-adjusted-age_smoothed-k9-1663_compiled-results.csv", row.names = 1) %>%
  mutate(sex= factor(sex, levels=c("F","M")),
         group= factor(group, levels=c("NTC.MDD","NTC.BPD","MDD.BPD")),
         dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased", "increased")))

restr.results_sm <- read.csv("processed-data/07_dx_DE/layer-restricted-age_smoothed-k9-1663_compiled-results.csv", row.names=1) %>%
  mutate(sex= factor(sex, levels=c("F","M")),
         group= factor(group, levels=c("NTC.MDD","NTC.BPD","MDD.BPD")),
         smoothed=factor(cluster, levels=c("L1","L2","L3.4","L5","L6","WM")),
         dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased", "increased")))


## Seurat label transfer
adj.results_se = read.csv("processed-data/07_dx_DE/layer-adjusted-age_seurat-pc30-no-lowUMI_compiled-results.csv", row.names=1) %>%
  mutate(sex= factor(sex, levels=c("F","M")),
         group= factor(group, levels=c("NTC.MDD","NTC.BPD","MDD.BPD")),
         dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased", "increased")))

restr.results_se <- read.csv("processed-data/07_dx_DE/layer-restricted-age_seurat-pc30-no-lowUMI_compiled-results.csv", row.names=1) %>%
  mutate(sex= factor(sex, levels=c("F","M")),
         group= factor(group, levels=c("NTC.MDD","NTC.BPD","MDD.BPD")),
         seurat_label_f=factor(cluster, levels=c("Micro.Vasc","Astro","L2.3","L4","L5","L6","Oligo","Inhb"),
                               labels=c("M.V","Astro","L2.3","L4","L5","L6","Oligo","Inhb")),
         dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased", "increased")))

# get lists of overlaps
checkList = list("NTC.MDD_F_dn"=c("group"="NTC.MDD","sex"="F","dir"="decreased"),
                 "NTC.MDD_F_up"=c("group"="NTC.MDD","sex"="F","dir"="increased"),
                 "NTC.MDD_M_dn"=c("group"="NTC.MDD","sex"="M","dir"="decreased"),
                 "NTC.MDD_M_up"=c("group"="NTC.MDD","sex"="M","dir"="increased"),
                 "NTC.BPD_F_dn"=c("group"="NTC.BPD","sex"="F","dir"="decreased"),
                 "NTC.BPD_F_up"=c("group"="NTC.BPD","sex"="F","dir"="increased"),
                 "NTC.BPD_M_dn"=c("group"="NTC.BPD","sex"="M","dir"="decreased"),
                 "NTC.BPD_M_up"=c("group"="NTC.BPD","sex"="M","dir"="increased"),
                 "MDD.BPD_F_dn"=c("group"="MDD.BPD","sex"="F","dir"="decreased"),
                 "MDD.BPD_F_up"=c("group"="MDD.BPD","sex"="F","dir"="increased"),
                 "MDD.BPD_M_dn"=c("group"="MDD.BPD","sex"="M","dir"="decreased"),
                 "MDD.BPD_M_up"=c("group"="MDD.BPD","sex"="M","dir"="increased")
)

restr_sm_tmp = bind_rows(filter(restr.results_sm, cluster!="WM", adj.P.Val<.05),
                         filter(restr.results_sm, cluster=="WM", sex=="M", adj.P.Val<.05),
                         filter(restr.results_sm, cluster=="WM", sex=="F", adj.P.Val<.01))
restr_se_tmp = bind_rows(filter(restr.results_se, cluster!="Oligo", adj.P.Val<.05),
                         filter(restr.results_se, cluster=="Oligo", sex=="M", adj.P.Val<.05),
                         filter(restr.results_se, cluster=="Oligo", sex=="F", adj.P.Val<.01))
checkList2 <- lapply(checkList, function(x) {
  outList = list()
  outList$adj_sm = filter(adj.results_sm, sex==x[["sex"]], group==x[["group"]], dir==x[["dir"]],
                          adj.P.Val<.05, abs(logFC)>.2)$gene_name
  outList$restr_sm = filter(restr_sm_tmp, sex==x[["sex"]], group==x[["group"]], dir==x[["dir"]],
                            abs(logFC)>.2)$gene_name
  outList$adj_se = filter(adj.results_se, sex==x[["sex"]], group==x[["group"]], dir==x[["dir"]],
                          adj.P.Val<.05, abs(logFC)>.2)$gene_name
  outList$restr_se = filter(restr_se_tmp, sex==x[["sex"]], group==x[["group"]], dir==x[["dir"]],
                            abs(logFC)>.2)$gene_name
  return(outList)
})

lapply(checkList2, function(x) sapply(x, length))


checkList3 <- lapply(checkList2, function(x) {
  outList = list()
  outList$sm_both = intersect(x$adj_sm, x$restr_sm)
  outList$se_both = intersect(x$adj_se, x$restr_se)
  return(outList)
})

lapply(checkList3, function(x) sapply(x, length))
#really interesting, all about equal
### MDD.BPD F up quite small but thats ok


checkList4 <- lapply(checkList3, function(x) intersect(x[[1]], x[[2]]))

sapply(checkList4, length)
#NTC.MDD_F_dn NTC.MDD_F_up NTC.MDD_M_dn NTC.MDD_M_up NTC.BPD_F_dn NTC.BPD_F_up 
#          29          102            0            0           51           22 
#NTC.BPD_M_dn NTC.BPD_M_up MDD.BPD_F_dn MDD.BPD_F_up MDD.BPD_M_dn MDD.BPD_M_up 
#.         31           33           17            3           33           45 


# save overlaps list
saveList <- list("sig_genes"= checkList2,
	"LA.LR_overlap"= checkList3,
	"LA.LR_both.annotations"= checkList4
)
saveRDS(saveList, "processed-data/07_dx_DE/LA-LR-overlap_lists.rds")
cat("\n\nOverlaps list saved to: processed-data/07_dx_DE/LA-LR-overlap_lists.rds\n")


# format for volcano plots
nlist = names(checkList2)
names(nlist) = nlist

volcano.df = do.call(rbind, lapply(nlist, function(x) {
  tmp = unlist(strsplit(x, "_"))
  target_group = tmp[[1]]
  target_sex = tmp[[2]]
  target_dir = tmp[[3]]
  target_dir2 = ifelse(target_dir=="dn", "decreased", "increased")
  
  ## precast smoothed
  la.ns_sm = filter(adj.results_sm, group == target_group, sex == target_sex, adj.P.Val>.05 | abs(logFC)<.02)
  la.filt_sm = filter(adj.results_sm, group == target_group, sex == target_sex) %>%
    mutate(filter_dir= target_dir2,
           is_sig_dir= gene_name %in% checkList2[[x]]$adj_sm, #grey, otherwise, black
           is_restr_sig = gene_name %in% checkList3[[x]]$sm_both, #light blue
           is_sig_both = gene_name %in% checkList4[[x]],
           point_color = factor(paste(is_sig_dir, is_restr_sig), levels=c("FALSE FALSE","TRUE FALSE","TRUE TRUE"),
                                labels=c("L-A NS", paste0("L-A only (", target_dir, ")"),
                                         paste0("L-A & L-R (", target_dir, ")"))),
           #both_color = factor(paste(is_sig_dir, is_restr_sig, is_sig_both), 
           #                     levels=c("FALSE FALSE FALSE","TRUE FALSE FALSE","TRUE TRUE FALSE", "TRUE TRUE TRUE"),
           #                     labels=c("L-A NS", "L-A only",
           #                              paste0("L-A & L-R (", target_dir,")"),
           #                              paste0("L-A & L-R (", target_dir,")\n(both annotations)"))),
           annotation="sm"
    )
  
  ## seurat labels
  la.filt_se = filter(adj.results_se, group == target_group, sex == target_sex) %>%
    mutate(filter_dir= target_dir2,
           is_sig_dir= gene_name %in% checkList2[[x]]$adj_se, #grey, otherwise, black
           is_restr_sig = gene_name %in% checkList3[[x]]$se_both, #light blue
           is_sig_both = gene_name %in% checkList4[[x]],
           point_color = factor(paste(is_sig_dir, is_restr_sig), levels=c("FALSE FALSE","TRUE FALSE","TRUE TRUE"),
                                labels=c("L-A NS", paste0("L-A only (", target_dir, ")"), 
                                         paste0("L-A & L-R (", target_dir, ")"))),
           #both_color = factor(paste(is_sig_dir, is_restr_sig, is_sig_both), 
           #                     levels=c("FALSE FALSE FALSE","TRUE FALSE FALSE","TRUE TRUE FALSE", "TRUE TRUE TRUE"),
           #                     labels=c("L-A NS", "L-A sig.",
           #                              paste0("L-A & L-R (", target_dir,")"),
           #                              paste0("L-A & L-R (", target_dir,")\n(both annotations)"))),
           annotation="se"
    )
  bind_rows(la.filt_sm, la.filt_se)
})
) %>% filter(dir==filter_dir)
rownames(volcano.df) <- NULL

table(volcano.df$point_color)


## make volcano plots
point_pal = c("L-A NS"="grey","L-A onl"="black","L-A only (dn)"="skyblue", "L-A only (up)"="tomato",
              "L-A & L-R (dn)"="skyblue4", "L-A & L-R (up)"="red3")


volcano.df_sm = filter(volcano.df, annotation=="sm")
volcano.df_se = filter(volcano.df, annotation=="se")

volcanoList <- list()
volcanoList$sm_single <- ggplot(volcano.df_sm, aes(x=logFC, y=-log10(adj.P.Val), color=point_color))+
  geom_point(shape=1)+scale_color_manual("", values=point_pal)+
  xlim(-max(volcano.df_sm$logFC),max(volcano.df_sm$logFC))+
  facet_grid(cols=vars(sex), rows=vars(group))+
  labs(title="PRECAST smoothed L-A results (age, PC3 covar)")+
  theme_bw()#+theme(legend.position="bottom")

volcanoList$se_single <- ggplot(volcano.df_se, aes(x=logFC, y=-log10(adj.P.Val), color=point_color))+
  geom_point(shape=1)+scale_color_manual("", values=point_pal)+
  xlim(-max(volcano.df_se$logFC),max(volcano.df_se$logFC))+
  facet_grid(cols=vars(sex), rows=vars(group))+
  labs(title="Seurat labels L-A results (age, PC3 covar)")+
  theme_bw()#+theme(legend.position="bottom")



volcanoList$sm_both <- ggplot(volcano.df_sm, aes(x=logFC, y=-log10(adj.P.Val), color=point_color, shape=is_sig_both))+
  geom_point()+scale_color_manual("", values=point_pal)+
  scale_shape_manual("Sig. both\nannotations", values=c(1,16))+
  xlim(-max(volcano.df_sm$logFC),max(volcano.df_sm$logFC))+
  facet_grid(cols=vars(sex), rows=vars(group))+
  labs(title="PRECAST smoothed L-A results (age, PC3 covar)")+
  theme_bw()

volcanoList$se_both <- ggplot(volcano.df_se, aes(x=logFC, y=-log10(adj.P.Val), color=point_color, shape=is_sig_both))+
  geom_point()+scale_color_manual("", values=point_pal)+
  scale_shape_manual("Sig. both\nannotations", values=c(1,16))+
  xlim(-max(volcano.df_se$logFC),max(volcano.df_se$logFC))+
  facet_grid(cols=vars(sex), rows=vars(group))+
  labs(title="Seurat labels L-A results (age, PC3 covar)")+
  theme_bw()



# format for bar plotting

ilist = lapply(nlist, function(x) {
  as.data.frame(list("coef_dir"=rep(x, 6),
                     "annotation"=c("sm","sm","sm","se","se","se"),
                     "model"=c("L-A only","L-R only","L-A & L-R",
                               "L-A only","L-R only","L-A & L-R"),
                     "set_size"=c(length(unique(setdiff(checkList2[[x]]$adj_sm, checkList3[[x]]$sm_both))),
                                  length(unique(setdiff(checkList2[[x]]$restr_sm, checkList3[[x]]$sm_both))),
                                  length(unique(checkList3[[x]]$sm_both)),
                                  length(unique(setdiff(checkList2[[x]]$adj_se, checkList3[[x]]$se_both))),
                                  length(unique(setdiff(checkList2[[x]]$restr_se, checkList3[[x]]$se_both))),
                                  length(unique(checkList3[[x]]$se_both))
                     )
  ))
})

df1 = do.call(rbind, ilist)
rownames(df1) <- NULL

df2 = as.data.frame(do.call(rbind, strsplit(df1$coef_dir, split="_")))
colnames(df2) = c("group","sex","dir")

df3 = mutate(cbind(df1, df2), coef=paste(group, sex))
df3$coef = factor(df3$coef, levels=c("NTC.MDD F","NTC.MDD M","NTC.BPD F","NTC.BPD M","MDD.BPD F","MDD.BPD M"))
df3$annotation = factor(df3$annotation, levels=c("sm","se"), labels=c("PRECAST (smoothed)", "Seurat labels"))




## make bar plots

df4 = filter(df3, model %in% c("L-A only", "L-A & L-R"))
df4.summ = group_by(df4, coef_dir, annotation) %>% mutate(total_set_size=sum(set_size)) %>%
  ungroup() %>%
  mutate(prop_set_size=set_size/total_set_size,
         mod_prop1 = ifelse(model=="L-A & L-R", .1, .9),
         mod_prop2= ifelse(prop_set_size==0, NA, mod_prop1))


fill_pal = c("L-A only (dn)"="skyblue","L-A only (up)"="tomato",
             "L-A & L-R (dn)"="skyblue4", "L-A & L-R (up)"="red3")

df4$fill_color= factor(paste(df4$model, factor(df4$dir, levels=c("dn","up"), labels=c("(dn)","(up)"))),
                       levels=names(fill_pal))



barList <- list()

df4_sm = filter(df4, annotation=="PRECAST (smoothed)")
df4.summ_sm = filter(df4.summ, annotation=="PRECAST (smoothed)")

barList$sm_la <- ggplot(df4_sm, aes(y=factor(coef, levels=rev(levels(coef))), fill=fill_color))+
  geom_bar(data=filter(df4_sm, dir=="dn"), aes(x= -set_size), stat="identity", position="fill", 
           color="black", linewidth=.3, width=.7)+
  geom_text(data=filter(df4.summ_sm, dir=="dn"), aes(x= -mod_prop2, label= set_size, fill=NULL))+
  geom_bar(data=filter(df4_sm, dir=="up"), aes(x= set_size), stat="identity", position="fill", 
           color="black", linewidth=.3, width=.7)+
  geom_text(data=filter(df4.summ_sm, dir=="up"), aes(x= mod_prop2, label= set_size, fill=NULL))+
  scale_fill_manual(values=fill_pal)+
  labs(title="PRECAST smoothed L-A results (age, PC3 covar)", x="prop. of sig. genes", fill="", y="")+
  theme_minimal()+theme(axis.text.y=element_text(size=12))

df4_se = filter(df4, annotation=="Seurat labels")
df4.summ_se = filter(df4.summ, annotation=="Seurat labels")

barList$se_la <- ggplot(df4_se, aes(y=factor(coef, levels=rev(levels(coef))), fill=fill_color))+
  geom_bar(data=filter(df4_se, dir=="dn"), aes(x= -set_size), stat="identity", position="fill", 
           color="black", linewidth=.3, width=.7)+
  geom_text(data=filter(df4.summ_se, dir=="dn"), aes(x= -mod_prop2, label= set_size, fill=NULL))+
  geom_bar(data=filter(df4_se, dir=="up"), aes(x= set_size), stat="identity", position="fill", 
           color="black", linewidth=.3, width=.7)+
  geom_text(data=filter(df4.summ_se, dir=="up"), aes(x= mod_prop2, label= set_size, fill=NULL))+
  scale_fill_manual(values=fill_pal)+
  labs(title="Seurat labels L-A results (age, PC3 covar)", x="prop. of sig. genes", fill="", y="")+
  theme_minimal()+theme(axis.text.y=element_text(size=12))


# table count genes

mytheme <- ttheme_default(
  core = list(fg_params=list(cex = .6)),
  colhead = list(fg_params=list(cex = .6)))

grobList <- lapply(nlist, function(x) {
  tmp = unlist(strsplit(x, "_"))
  target_group = tmp[[1]]
  target_sex = tmp[[2]]
  target_dir = tmp[[3]]
  target_dir2 = ifelse(target_dir=="dn", "decreased", "increased")
  
  t1 = as.data.frame(list("set"=c("L-A","L-A & L-R"),
                          "sm"=c(length(checkList2[[x]]$adj_sm), length(checkList3[[x]]$sm_both)),
                          "se"=c(length(checkList2[[x]]$adj_se), length(checkList3[[x]]$se_both)))
  )
  gt1 = tableGrob(t1, rows = NULL, cols=c(paste("Sig.", target_dir2, sep="\n"),
                                          "PRECAST\n(smoothed)","Seurat labels"),
                  theme = mytheme)
  
  if(target_dir=="dn") fill_color = "skyblue"
  if(target_dir=="up") fill_color = "tomato"
  t2 = as.data.frame(list("set"="L-A & L-R\n(both annotations)",
                          "sm"=length(checkList4[[x]])))
  tt2 <- ttheme_default(core=list(bg_params = list(fill=c(fill_color,fill_color)),
                                  fg_params=list(cex = .6)),
                        colhead = list(fg_params=list(cex = .6))
                        )
  gt2 = tableGrob(t2, rows = NULL, cols=NULL, theme=tt2)
  
  gt3 = gtable_combine(gt1, gt2, along=2, join="left")
  gt3$layout[c(20,22),"r"] = 3
  
  return(gt3)
})


## combine directions from same group x sex

tmp.mtx = matrix(1:length(nlist), ncol=2, byrow=T)

grobList2 = lapply(1:nrow(tmp.mtx), function(x) {
  n1 = names(grobList)[[tmp.mtx[x,1]]]
  n2 = names(grobList)[[tmp.mtx[x,2]]]
  tmp = unlist(strsplit(n1, "_"))
  target_group = tmp[[1]]
  target_sex = tmp[[2]]
  
  t1 <- gtable_combine(grobList[[n2]], grobList[[n1]], along=2)
  title <- textGrob(paste(target_group, target_sex),gp=gpar(fontsize=20))
  padding <- unit(5,"mm")
  
  table <- gtable_add_rows(
    t1, 
    heights = grobHeight(title) + padding,
    pos = 0)
  table <- gtable_add_grob(
    table, 
    title, 
    1, 1, 1, ncol(table))
  
  return(table)
})



# compile pdf
pdf(file="plots/07_dx_DE/LA-LR-overlap_summary.pdf", height=11, width=8)
#p1
grid.arrange(rasterize(volcanoList$sm_single+theme(plot.margin = ggplot2::margin(1,1,1,1, unit="cm")), dpi=200), 
             barList$sm_la+theme(aspect.ratio=1), 
             ncol=1, top="PRECAST (smoothed)")
#p2
grid.arrange(rasterize(volcanoList$se_single+theme(plot.margin = ggplot2::margin(1,1,1,1, unit="cm")), dpi=200), 
             barList$se_la+theme(aspect.ratio=1), ncol=1, top="Seurat labels")
#p3
grid.arrange(grobList2[[1]], grobList2[[2]],
            grobList2[[3]], grobList2[[4]],
            grobList2[[5]], grobList2[[6]], 
            layout_matrix=matrix(1:6, ncol=2, byrow=T)
)
#p4
grid.arrange(rasterize(volcanoList$sm_both+theme(plot.margin = ggplot2::margin(1,1,1,1, unit="cm")), dpi=200), 
             rasterize(volcanoList$se_both+theme(plot.margin = ggplot2::margin(1,1,1,1, unit="cm")), dpi=200), ncol=1, top="L-A & L-R in both annotations")
dev.off()
cat("\n\nSaved L-A and L-R overlap strategy summary to: plots/07_dx_DE/LA-LR-overlap_summary.pdf\n")


cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
