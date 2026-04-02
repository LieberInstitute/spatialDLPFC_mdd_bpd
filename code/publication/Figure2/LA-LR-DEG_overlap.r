setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
})

cpList <- readRDS("plots/colorPalettes.rds")

comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")
names(comparisons) <- comparisons
comparisons2 = comparisons[c(1,3,5,2,4,6)]

la.sm = read.csv("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv")
lr.sm = read.csv("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv")

is_sig = function(x) x!="NS"

comb_names = function(x, y) paste(sort(c(x,y)), collapse=" - ")

laList = lapply(comparisons, function(i) {
  tmp.la = la.sm[,c("gene_id", "gene_name", paste0(i, "_coef"), paste0(i, "_ttest"))]
  colnames(tmp.la)[3:4] = c("dir_LA","sig_LA")
  filter(tmp.la, sig_LA!="NS") %>% mutate(dir_LA=sign(dir_LA), sex.group=i)
})

cat("\nAny L-A gene sig. in different direction in different dx*sex groups?\n")
do.call(rbind, laList) %>% filter(sex.group %in% c("F_NTC.MDD", "M_NTC.MDD", "F_NTC.BPD", "M_NTC.BPD")) %>%
  group_by(gene_id, gene_name, dir_LA) %>% summarise(sig.groups=paste(sex.group, collapse="/")) %>%
  group_by(gene_id, gene_name) %>% add_tally(name="n_dir") %>%
  filter(n_dir>1)
cat("\n")
# there are two genes that are sig in different directions
# CIRBP decreased F_NTC.MDD and increased M_NTC.MDD and M_NTC.BPD
filter(la.sm, gene_name=="CIRBP") # less than .2 logFC all around
cat("\n")
# RBM3 decreased F_NTC.MDD and F_NTC.BPD and increased M_NTC.BPD
filter(la.sm, gene_name=="RBM3") # greater than .2 logFC all around
cat("\n")

lrList = lapply(comparisons, function(i) {
  cnames = grep(paste0(i,"_ttest"), colnames(lr.sm), value=T)
  clusnames = gsub(paste0("_",i,"_ttest"), "", cnames)
  
  sig.dir = sign(lr.sm[,grep(paste0(i,"_coef"), colnames(lr.sm))])
  colnames(sig.dir) = clusnames
  sig.status = lr.sm[,cnames]!="NS"
  colnames(sig.status) = clusnames
  
  sig.dir[!sig.status] = 0
  
  tmp = cbind(select(lr.sm, gene_id, gene_name, "n_ttest_sig"=paste0("n_ttest_sig_",i)),
              sig.dir) %>%
    filter(n_ttest_sig>0) %>% mutate(sex.group=i)
  
  return(tmp)
})


cat("\nAny L-R gene significant in different direction in different domains?\n")
sapply(lrList, function(x) {
 check_condition = apply(x[,4:9], MARGIN=1, function(x) length(sign(unique(x[x!=0])))>1) 
 sum(check_condition)
})

cat("\nAny L-R gene sig. in different direction in different dx*sex groups?\n")
m1 = do.call(rbind, lapply(lrList, function(x) {
  data.frame("sex.group"=x$sex.group, "gene_name"=x$gene_name, 
             "dir"=apply(x[,4:9], MARGIN=1, function(x) sign(unique(x[x!=0]))))
})) %>%
  tidyr::pivot_wider(names_from="sex.group", values_from="dir")
check_condition = apply(m1[,2:5], MARGIN=1, function(x) length(sign(unique(x[!is.na(x)])))>1) 
m1[check_condition,]
cat("\n")
do.call(rbind, lrList) %>% filter(gene_name=="RBM3")
cat("\n")


overlapList <- lapply(comparisons, function(i) {
  cnames = grep(paste0(i,"_ttest"), colnames(lr.sm), value=T)
  clusnames = gsub(paste0("_",i,"_ttest"), "", cnames)
  
  tmp.lr = lrList[[i]]; tmp.la = laList[[i]];
  tmp.lr = tidyr::pivot_longer(tmp.lr, all_of(clusnames), names_to="cluster", values_to="dir_LR") %>%
    filter(dir_LR!=0)
  
  tmp.comb = left_join(tmp.lr, tmp.la[,1:3]) %>%
    mutate(dir_LA=ifelse(is.na(dir_LA), 0, dir_LA),
           comparison=i)
  
  out1 = group_by(tmp.comb, comparison, cluster, dir_LA, dir_LR) %>% tally(name="n_DEGs") %>%
    ungroup()
  
  if(length(setdiff(clusnames, out1$cluster))>0) {
    out1 = add_row(out1, comparison=i, cluster=setdiff(clusnames, out1$cluster),
                   dir_LA=0, dir_LR=1, n_DEGs=0) %>%
      add_row(comparison=i, cluster=setdiff(clusnames, out1$cluster),
              dir_LA=1, dir_LR=1, n_DEGs=0)
  }
  
  if(length(setdiff(tmp.la$gene_id, tmp.comb$gene_id))>0) {
    la.only = filter(tmp.la, gene_id %in% setdiff(tmp.la$gene_id, tmp.comb$gene_id))
    out2 = group_by(la.only, dir_LA) %>% tally()
    out1 = add_row(out1, comparison=i, cluster="L-A", 
                   dir_LA=out2$dir_LA, dir_LR=0, n_DEGs=out2$n)
  }
  
  return(mutate(out1, cluster=factor(cluster, levels=c("L-A",clusnames))))
})

## in order for the x axis to be appropriately labeled in the plot, I need to 
## make an entry for all possible options for F_NTC.MDD (the first group plotted)
### if I make a different plot for each comparison (rather than faceting by comparison), 
### i will need to do this for all comparisons
#names(overlapList)
#filter(overlapList[[1]], cluster!="L-A") %>% group_by(cluster) %>% add_tally() %>% filter(n<4)
overlapList[[1]] = add_row(overlapList[[1]], comparison="F_NTC.MDD", cluster=c("L5","L5","L6","L6"),
         dir_LA=c(0,0,0,0), dir_LR=c(-1,1,-1,1), n_DEGs=c(0,0,0,0))

overlap.df = do.call(rbind, overlapList) %>% 
  mutate(comparison=factor(comparison, levels=names(comparisons2)),
         cluster= factor(cluster, levels=c("L-A", names(cpList$smoothed.bright))))

cat("\nAny gene sig. for both L-A and L-R but in different directions?\n")
filter(overlap.df, dir_LA!=0 & dir_LR!=0, sign(dir_LA)!=sign(dir_LR))
cat("\n")


# format for bar plot
overlap.df = mutate(overlap.df, n_models= ifelse(dir_LA!=0 & dir_LR!=0, "two models", "one model"),
       dir_overall= ifelse(dir_LR<0 | dir_LA<0, "decreased", "increased"),
       n_DEGs_dir = ifelse(dir_overall=="decreased", -n_DEGs, n_DEGs))

overlap.df$facets = factor(paste(overlap.df$cluster, overlap.df$n_models), levels=c("L-A one model",
                                                paste(names(cpList$smoothed.bright), "two models"),
                                                paste(names(cpList$smoothed.bright), "one model")),
                     labels=c("L-A\nonly", rep("L-A & L-R\noverlap", length(cpList$smoothed.bright)),
                              rep("L-R only", length(cpList$smoothed.bright))))

p1 <- ggplot(overlap.df, aes(x=cluster, y=n_DEGs_dir, fill=cluster))+
  geom_bar(data=filter(overlap.df, dir_overall=="decreased"), aes(group=n_models), 
           stat="identity", position="stack", color="black", linewidth=.3, width=.7)+
  geom_bar(data=filter(overlap.df, dir_overall=="increased"), aes(group=n_models), 
           stat="identity", position="stack", color="black", linewidth=.3, width=.7)+
  geom_hline(aes(yintercept=0), linewidth=1)+
  scale_y_continuous(breaks=c(-40,-20,0,20,40), limits=c(-50,50))+
  scale_fill_manual(values=c("L-A", cpList$smoothed.bright), guide="none")+
  facet_grid(cols=vars(facets), rows=vars(comparison), scales="free_x", space="free")+
  labs(y="# DEGs (directional)", x="domain")+
  theme_bw()+theme(strip.background = element_rect(fill="white", color=NA),
                   panel.grid.major.x=element_blank(), panel.grid.minor=element_blank())


# format for L-R between domain overlap heatmap
compList = lapply(comparisons, function(i) {
  cnames = grep(paste0(i,"_ttest"), colnames(lr.sm), value=T)
  clusnames = gsub(paste0("_",i,"_ttest"), "", cnames)
  
  tmp = cbind(select(lr.sm, gene_id, gene_name, "n_ttest_sig"=paste0("n_ttest_sig_",i)),
              lr.sm[,cnames]) %>%
    filter(n_ttest_sig>0) %>%
    mutate_at(cnames, is_sig)
  colnames(tmp)[4:ncol(tmp)] <- clusnames
  
  plot.df = do.call(rbind, lapply(clusnames, function(j) {
    t1 = tmp[tmp[[j]],]
    if(nrow(t1)==0) {
      out1 = data.frame("query"=clusnames, n_overlap=0, reference=j, n_ref_total=0)
    } else {
      t2 = tidyr::pivot_longer(t1, setdiff(clusnames, j), names_to="query", values_to="is_sig") %>%
        group_by(query) %>% summarise(n_overlap=sum(is_sig))
      lr.unique = t1$gene_name[rowSums(t1[,setdiff(clusnames, j)])==0]
      out1 = add_row(t2, query=j, n_overlap=length(lr.unique)) %>%
        mutate(reference=j, n_ref_total=nrow(t1))
    }
    out1$clus_pair = sapply(seq(nrow(out1)), function(z) comb_names(out1$query[[z]], out1$reference[[z]]))
    return(out1)
  })) %>% mutate(reference=factor(reference, levels=rev(clusnames)),
                 query=factor(query, levels=clusnames), 
                 comparison=i)
  
  return(plot.df)
})


comp.df = do.call(rbind, compList) %>% 
  mutate(fill_bin = cut(n_overlap, breaks=c(0,1,4,7,10,13,60), include.lowest = T, right=F,
                        labels=c("0","1-3","4-6","7-9","9-12",">12")),
         comparison=factor(comparison, levels=names(comparisons2)))

p2 <- ggplot(mutate(comp.df, facets="col\n1"), aes(x=query, y=reference, fill=fill_bin))+
  geom_tile(color="black", linewidth=.3)+
  scale_fill_manual(values=c("white",RColorBrewer::brewer.pal(n=5,"Greys")[2:5],"black"))+
  facet_grid(rows=vars(comparison), cols=vars(facets))+
  labs(x="domain", fill="# lrDEGs\noverlap")+
  theme_minimal()+theme(legend.position="none", axis.title.y=element_blank(),
                        strip.background = element_rect(fill="white", color=NA),
                        panel.grid.major=element_blank(), panel.grid.minor=element_blank())

# separate plot for legend
plegend <- ggplot(comp.df, aes(x=query, y=reference, fill=fill_bin))+
  geom_tile(color="black", linewidth=.3)+
  scale_fill_manual(values=c("white",RColorBrewer::brewer.pal(n=5,"Greys")[2:5],"black"))+
  facet_wrap(vars(comparison), ncol=2)+
  labs(x="domain", y="domain", fill="# lrDEGs\noverlap")+
  theme_minimal()

pdf(file="plots/publication/Figure2/LA-LR-DEG_overlap.pdf", height=6, width=6)
gridExtra::grid.arrange(p1, p2, layout_matrix=matrix(c(1,1,1,2), ncol=4))
plegend
dev.off()

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
