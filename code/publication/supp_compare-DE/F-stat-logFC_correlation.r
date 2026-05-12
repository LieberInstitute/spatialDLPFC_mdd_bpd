setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
	library(ggrastr)
})

resultsList <- c("smoothed-k9-1663","seurat-pc30")
comp_names <- list("smoothed-k9-1663"= c("L1","L2","L3.4","L5","L6","WM"),
                   "seurat-pc30"=  c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))

comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")
names(comparisons) <- comparisons
comparisons2 = comparisons[c(1,3,5,2,4,6)]


## L-A
# the layer-adjusted model is asking fundamentally the same question between the two annotations
sm.la = read.csv("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_smoothed-k9-1663_rev-gene-input_F-test.csv")
se.la = read.csv("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_seurat-pc30_rev-gene-input_F-test.csv")

cat("\n*** Layer-adjusted F-test overlap ***")
cat("\nBased on all genes tested for each annotation...")
sm.la.sig = filter(sm.la, adj.P.Val<.05)$gene_name
cat("\n# of PRECAST sig:",length(sm.la.sig))
se.la.sig = filter(se.la, adj.P.Val<.05)$gene_name
cat("\n# of Seurat sig:", length(se.la.sig))

cat("\n# of overlap sig:", length(intersect(sm.la.sig, se.la.sig)))
cat("\n# of only PRECAST sig:", length(setdiff(sm.la.sig, se.la.sig)))
cat("\n# of only Seurat sig:", length(setdiff(se.la.sig, sm.la.sig)))

both.genes = intersect(sm.la$gene_id, se.la$gene_id)
cat("\n\nBased on genes tested in both annotations...")
cat("\nBoth genes:", length(both.genes)) #13015
rownames(sm.la) <- sm.la$gene_id
rownames(se.la) <- se.la$gene_id


sm.la.sig2 = filter(sm.la[both.genes,], adj.P.Val<.05)$gene_name
cat("\n# of PRECAST sig:",length(sm.la.sig2))
se.la.sig2 = filter(se.la[both.genes,], adj.P.Val<.05)$gene_name
cat("\n# of Seurat sig:",length(se.la.sig2))

cat("\n# of overlap sig:", length(intersect(sm.la.sig2, se.la.sig2)))
cat("\n# of only PRECAST sig:", length(setdiff(sm.la.sig2, se.la.sig2)))
cat("\n# of only Seurat sig:", length(setdiff(se.la.sig2, sm.la.sig2)))


colnames(sm.la)[c(8,10)] = c("F_sm","adj.P.Val_sm")
colnames(se.la)[c(8,10)] = c("F_se","adj.P.Val_se")

la.f = merge(sm.la[both.genes,c("gene_id","gene_name","F_sm","adj.P.Val_sm")],
      se.la[both.genes,c("gene_id","gene_name","F_se","adj.P.Val_se")])

la.f.corr = round(cor(la.f$F_sm, la.f$F_se), 4)
p1 <- ggplot(la.f, aes(x=F_sm, y=F_se))+
  rasterize(geom_point(size=.1), dpi=300)+
  geom_vline(aes(xintercept=min(filter(la.f, adj.P.Val_sm<.05)$F_sm)), lty=2, color="red")+
  geom_hline(aes(yintercept=min(filter(la.f, adj.P.Val_se<.05)$F_se)), lty=2, color="red")+
  coord_cartesian(ylim=c(0,25), xlim=c(0,25))+
  labs(title=paste0("rho= ", la.f.corr), 
       x="PRECAST domains", y="Seurat labels")+
  theme_bw()+theme(panel.grid.minor=element_blank(), aspect.ratio=1,
	text=element_text(size=8))


df1 <- do.call(rbind, lapply(comparisons2, function(x) {
  tmp = data.frame(sex.group=x,"logFC_sm"=sm.la[both.genes,x], "logFC_se"=se.la[both.genes,x])
  cor.tmp = round(cor(tmp$logFC_sm, tmp$logFC_se), 4)
  tmp$rho = cor.tmp
  return(tmp)
})) %>% mutate(sex.group=factor(sex.group, levels=comparisons2))

facet_labels = distinct(df1, sex.group, rho)
facet_labels = paste(facet_labels$sex.group, "rho=", round(facet_labels$rho,2))
names(facet_labels) = names(comparisons2)
p2 <- ggplot(df1, aes(x=logFC_sm, y=logFC_se))+
  rasterize(geom_point(size=.1), dpi=300)+
  facet_wrap(vars(sex.group), ncol=6, labeller=as_labeller(facet_labels))+
  ylim(-3,3)+xlim(-3,3)+
  labs(x="logFC PRECAST domains", y="logFC Seurat labels")+
  theme_bw()+theme(panel.grid.minor=element_blank(), text=element_text(size=8),
                   aspect.ratio=1,
                   strip.background = element_rect(fill="transparent", color=NA))

pdf(file="plots/publication/supp_compare-DE/F-stat-logFC_layer-adjusted_correlation.pdf", width=6.5, height=2.5)
grid.arrange(p1, p2, layout_matrix=matrix(c(1,2,2,2,2,2,2), ncol=7))
dev.off()


## L-R
# the layer-restricted model is asking fundamentally the different questions because the clusters are different
sm.lr = read.csv("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_smoothed-k9-1663_rev-gene-input_F-test_all.csv")
se.lr = read.csv("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_seurat-pc30_rev-gene-input_F-test_all.csv")


cat("\n\n\n*** Layer-restricted F-test overlap ***")
cat("\nBased on all genes tested for each annotation...")
sm.lr.sig = filter(sm.lr, adj.P.Val<.05)$gene_name
cat("\n# of PRECAST sig:",length(sm.lr.sig))
se.lr.sig = filter(se.lr, adj.P.Val<.05)$gene_name
cat("\n# of Seurat sig:", length(se.lr.sig))

cat("\n# of overlap sig:", length(intersect(sm.lr.sig, se.lr.sig)))
cat("\n# of only PRECAST sig:", length(setdiff(sm.lr.sig, se.lr.sig)))
cat("\n# of only Seurat sig:", length(setdiff(se.lr.sig, sm.lr.sig)))

both.genes = intersect(sm.lr$gene_id, se.lr$gene_id)
cat("\n\nBased on genes tested in both annotations...")
cat("\nBoth genes:", length(both.genes)) #13015
rownames(sm.lr) <- sm.lr$gene_id
rownames(se.lr) <- se.lr$gene_id


sm.lr.sig2 = filter(sm.lr[both.genes,], adj.P.Val<.05)$gene_name
cat("\n# of PRECAST sig:",length(sm.lr.sig2))
se.lr.sig2 = filter(se.lr[both.genes,], adj.P.Val<.05)$gene_name
cat("\n# of Seurat sig:",length(se.lr.sig2))

cat("\n# of overlap sig:", length(intersect(sm.lr.sig2, se.lr.sig2)))
cat("\n# of only PRECAST sig:", length(setdiff(sm.lr.sig2, se.lr.sig2)))
cat("\n# of only Seurat sig:", length(setdiff(se.lr.sig2, sm.lr.sig2)))
cat("\n\n")

colnames(sm.lr)[c(38,40)] = c("F_sm","adj.P.Val_sm")
colnames(se.lr)[c(50,52)] = c("F_se","adj.P.Val_se")

lr.f = merge(sm.lr[both.genes,c("gene_id","gene_name","F_sm","adj.P.Val_sm")],
             se.lr[both.genes,c("gene_id","gene_name","F_se","adj.P.Val_se")])

lr.f.corr = round(cor(lr.f$F_sm, lr.f$F_se), 4)
p3 <- ggplot(lr.f, aes(x=F_sm, y=F_se))+
  rasterize(geom_point(size=.1), dpi=300)+
  geom_vline(aes(xintercept=min(filter(lr.f, adj.P.Val_sm<.05)$F_sm)), lty=2, color="red")+
  geom_hline(aes(yintercept=min(filter(lr.f, adj.P.Val_se<.05)$F_se)), lty=2, color="red")+
  coord_cartesian(ylim=c(0,5.5), xlim=c(0,5.5))+
  labs(title=paste0("rho= ", lr.f.corr), 
       x="PRECAST domains", y="Seurat labels")+
  theme_bw()+theme(panel.grid.minor=element_blank(), aspect.ratio=1,
	text=element_text(size=8))


# logFC
df2 <- do.call(rbind, lapply(comparisons2, function(i) {
  sm1 = sm.lr[both.genes,grep(i, colnames(sm.lr))]
  colnames(sm1) <- paste0("sm_", colnames(sm1))
  se1 = se.lr[both.genes,grep(i, colnames(se.lr))]
  colnames(se1) <- paste0("se_", colnames(se1))
  stopifnot(identical(rownames(sm1), rownames(se1)))
  
  corr.m = cor(cbind(sm1, se1), method="pearson")
  corr.m2 = corr.m[grep("sm_", rownames(corr.m)), grep("se_", colnames(corr.m))]
  rownames(corr.m2) <- gsub(paste0("_",i), "", gsub("sm_","", rownames(corr.m2)))
  colnames(corr.m2) <- gsub(paste0("_",i), "", gsub("se_","", colnames(corr.m2)))
  
  cor.df = tidyr::pivot_longer(tibble::rownames_to_column(as.data.frame(corr.m2), var="PRECAST"), all_of(comp_names[[2]]),
                               names_to="clusters", values_to="pearson_r") %>%
    mutate(clusters=factor(clusters, levels=rev(comp_names[[2]]), labels=c("Olg","L6","L5","Inb","L4","L2.3","Ast","M.V")), 
           PRECAST=factor(PRECAST, levels=comp_names[[1]]),
           sex.group=i)
  return(cor.df)
})) %>% mutate(sex.group=factor(sex.group, levels=comparisons2))


p4 <- ggplot(df2, aes(x=PRECAST, y=clusters, fill=pearson_r))+
  geom_tile(color="grey50", linewidth=.1)+
  scale_fill_gradientn("Pearson\nrho", 
                       colors=colorRampPalette(RColorBrewer::brewer.pal(n=7, "RdBu")[7:1])(20),
                       limits=c(0, 1))+
  facet_wrap(vars(sex.group), ncol=6)+
  labs(x="logFC PRECAST domains", y="logFC Seurat labels",
       title="Layer-restricted logFC correlations")+
  theme_minimal()+theme(aspect.ratio=1, panel.grid = element_blank(),
                          text=element_text(size=8), plot.title=element_text(hjust=.5),
                           #axis.title.x=element_blank(), axis.title.y=element_blank(),
                          axis.ticks = element_line(color="grey50", linewidth=.3), axis.text.y=element_text(color="black"),
                          )

pdf(file="plots/publication/supp_compare-DE/F-stat-logFC_layer-restricted_correlation.pdf", width=6.5, height=2.5)
grid.arrange(p3, p4+theme(legend.position="none"), layout_matrix=matrix(c(1,2,2,2,2,2,2), ncol=7))
grid.arrange(p3, p4, layout_matrix=matrix(c(1,2,2,2,2,2,2), ncol=7))
dev.off()


cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
