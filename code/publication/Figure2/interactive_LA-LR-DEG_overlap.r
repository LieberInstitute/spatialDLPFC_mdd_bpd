setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  #library(SpatialExperiment)
  #library(scater)
  library(gridExtra)
})
set.seed(123)

cpList <- readRDS("plots/colorPalettes.rds")
source("code/06_pseudobulk/custom_functions.r")

results_set = "smoothed-k9-1663"
comp_names = c("L1","L2","L3dot4","L5","L6","WM")
names(comp_names) = c("L1","L2","L3.4","L5","L6","WM")
col.pal = cpList$smoothed.bright

comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")
names(comparisons) <- comparisons
comparisons2 = comparisons[c(1,3,5,2,4,6)]

la.degs = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))
lr.degs = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))

la.degs_t = filter(la.degs, n_ttest_sig>0)
lr.degs_t = filter(lr.degs, n_ttest_sig>0)


lr.and.la_t = intersect(la.degs_t$gene_name, lr.degs_t$gene_name)

la.not.lr_t = setdiff(la.degs_t$gene_name, lr.degs_t$gene_name)

lr.not.la_t = setdiff(lr.degs_t$gene_name, la.degs_t$gene_name)

la.unique = sapply(comparisons, function(i) {
  any.clus = lr.degs_t[[paste0("n_ttest_sig_",i)]]>0
  la.res = la.degs_t[[paste0(i,"_ttest")]]!="NS"
  c("L-A_L-R_overlap"=length(intersect(lr.degs_t$gene_name[any.clus], la.degs_t$gene_name[la.res])),
    "LA_only"=length(setdiff(la.degs_t$gene_name[la.res], lr.degs_t$gene_name[any.clus])),
    "LR_only"=length(setdiff(lr.degs_t$gene_name[any.clus], la.degs_t$gene_name[la.res])))
})


lr.unique = lapply(comparisons, function(i) {
  any.clus = lr.degs_t[[paste0("n_ttest_sig_",i)]]>0
  la.res = la.degs_t[[paste0(i,"_ttest")]]!="NS"
  list("LR_only"=setdiff(lr.degs_t$gene_name[any.clus], la.degs_t$gene_name[la.res]),
       "LR_LA_overlap"=intersect(lr.degs_t$gene_name[any.clus], la.degs_t$gene_name[la.res]))
})

lrt = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_",
                      results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
  mutate(dir=factor(sign(logFC), levels=c(1,-1), labels=c("Inc.","Dec.")), 
         cluster=factor(cluster, levels=names(comp_names)),
         sex.group = factor(paste(sex, group, sep="_"), levels=comparisons)) 

i="F_NTC.MDD"

c1 = cpList$smoothed.light
names(c1) = paste(names(c1), "L-A & L-R overlap")
c2 = cpList$smoothed.bright
names(c2) = paste(names(c2), "L-R only")
col.pal = c("L-A L-A only"="grey65", c1, c2)
grep("L-A & L-R overlap", names(col.pal))

c3 = cpList$smoothed.bright
names(c3) = paste(names(c3), "L-A & L-R overlap")
c4 = rep("black", length(cpList$smoothed.light))
names(c4) = paste(names(cpList$smoothed.light), "L-R only")
col.pal2 = c("L-A L-A only"="grey30", c3, c4)


plist <- lapply(comparisons2, function(i) {
  tmp = filter(lrt, sex.group==i, adj.P.Val<.05)
  
  tmp.lr.overlap = filter(tmp, gene_name %in% lr.unique[[i]]$LR_LA_overlap) %>%
    group_by(cluster, .drop=F) %>% tally() %>% mutate(overlap="L-A & L-R overlap")
  tmp.lr.only = filter(tmp, gene_name %in% lr.unique[[i]]$LR_only) %>% 
    group_by(cluster) %>% tally() %>% mutate(overlap="L-R only")
  tmp.la.only = data.frame(cluster="L-A", n=la.unique["LA_only",i], overlap="L-A only")
  
  tmp = bind_rows(tmp.lr.overlap, tmp.lr.only, tmp.la.only) %>% 
    mutate(cluster=factor(cluster, levels=c("L-A", levels(tmp.lr.overlap$cluster))),
           ocolor=paste(as.character(cluster), overlap),
           overlap=factor(overlap, levels=c("L-A only","L-A & L-R overlap","L-R only")))
  
  if(results_set=="smoothed-k9-1663") ymax=55
  if(results_set=="seurat-pc30") ymax=65
  if(i=="M_MDD.BPD") ymax=80
  
  p <- ggplot(tmp, aes(x=cluster, y=n, fill=ocolor, color=ocolor))+
    geom_bar(aes(group=overlap), stat="identity", position=position_dodge2(preserve="single"),
             linewidth=.3)+
    ylim(0,ymax)+
    scale_fill_manual(values=col.pal, guide="none")+
    scale_color_manual(values=col.pal2, guide="none")+
    #scale_color_manual(values=c("black","black","grey50"), guide="none")+
    labs(title=i, x="model", y="# DEGs")+
    theme_minimal()+theme(panel.border = element_rect(fill=NA, color="grey"),
                          plot.title=element_text(hjust=.5),
                          panel.grid.minor=element_blank(),
                          text=element_text(size=10))

  if(results_set=="seurat-pc30") p <- p+scale_x_discrete(labels=c("M.V","A","L2.3","L4","In","L5","L6","Olg"))+theme(axis.text.x=element_text(size=8))
  return(p)
})


l.df = mutate(data.frame(model=c("L-A", rep(c("L1","L2","L3.4","L5","L6","WM"), 2)),
           overlap=c("L-A only", rep("L-A & L-R overlap", 6), rep("L-R only", 6)),
           xpos=c("b",rep("a",6),rep("b",6))),
       ocolor= paste(model, overlap),
       model=factor(model, levels=rev(c("L-A", names(cpList$smoothed.light))))
)
lplot = ggplot(l.df, aes(x=xpos, y=model))+
  geom_count(aes(color=ocolor, fill=ocolor), shape=22, size=8)+
  geom_text(data=data.frame(xpos=c("b"), model=c("L1"), overlap=c("L-R only")), 
            aes(label=overlap), nudge_y=.5, nudge_x=-.2, hjust=0, vjust=.5, size=4)+
  geom_text(data=data.frame(xpos=c("a"), model=c("L1"), overlap=c("L-A & L-R overlap")), 
            aes(label=overlap), nudge_y=.5, nudge_x=.2, hjust=1, vjust=.5, size=4)+
  geom_text(data=data.frame(xpos=c("b",rep("a", 6)), model=c("L-A","L1","L2","L3.4","L5","L6","WM"), 
                            overlap=c("L-A only","L1","L2","L3.4","L5","L6","WM")), 
            aes(label=overlap), hjust=1, nudge_x=-.5, size=6)+
  scale_fill_manual(values=col.pal, guide="none")+
  scale_color_manual(values=col.pal2, guide="none")+
  scale_x_discrete(expand=expansion(add = c(2.5)))+
  scale_y_discrete(expand=expansion(add = c(1)))+
  theme_void()+theme(plot.margin=margin(.1,4,.1,4,"cm"), legend.position="none")

pdf(file="plots/publication/Figure2/LA-LR-DEG_overlap.pdf", width=6, height=5)
do.call(grid.arrange, c(plist, ncol=3))
lplot
dev.off()
