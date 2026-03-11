library(dplyr)
library(ggplot2)

resultsList <- c("smoothed-k9-1663","seurat-pc30")
comp_names <- list("smoothed-k9-1663"= c("L1","L2","L3.4","L5","L6","WM"),
                   "seurat-pc30"=  c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))

comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")
names(comparisons) <- comparisons

# all F test ----
## L-A ----
# the layer-adjusted model is asking fundamentally the same question between the two annotations
sm.la = read.csv("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_smoothed-k9-1663_rev-gene-input_F-test.csv")
se.la = read.csv("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_seurat-pc30_rev-gene-input_F-test.csv")
#cc.la = read.csv("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_custom-cluster_rev-gene-input_F-test.csv")

both.genes = intersect(sm.la$gene_id, se.la$gene_id)
length(both.genes) #13015
#all.genes = intersect(intersect(sm.la$gene_id, se.la$gene_id), cc.la$gene_id)
#length(all.genes) #12845
rownames(sm.la) <- sm.la$gene_id
rownames(se.la) <- se.la$gene_id
#rownames(cc.la) <- cc.la$gene_id

par(mfrow=c(1,1))
plot(sm.la[both.genes,"F"], se.la[both.genes,"F"], xlab="PRECAST F statistic", ylab="Seurat F statistic",
     main="Layer-adjusted", sub=paste("Pearson rho =",round(cor(sm.la[both.genes,"F"], se.la[both.genes,"F"], method="pearson"),3)))
cor(sm.la[both.genes,"F"], se.la[both.genes,"F"], method="pearson") #0.9545489

plot(cc.la[all.genes,"F"], sm.la[all.genes,"F"], xlab="Custom F statistic", ylab="PRECAST F statistic",
     main="Layer-adjusted", sub=paste("Pearson rho =",round(cor(cc.la[all.genes,"F"], sm.la[all.genes,"F"], method="pearson"),3)))
cor(cc.la[all.genes,"F"], sm.la[all.genes,"F"], method="pearson") #0.9402373

plot(cc.la[all.genes,"F"], se.la[all.genes,"F"], xlab="Custom F statistic", ylab="Seurat F statistic",
     main="Layer-adjusted", sub=paste("Pearson rho =",round(cor(cc.la[all.genes,"F"], se.la[all.genes,"F"], method="pearson"),3)))
cor(cc.la[all.genes,"F"], se.la[all.genes,"F"], method="pearson") #0.9771998

par(mfrow=c(3,2))
for(i in comparisons) {
  c1= cor(sm.la[both.genes,i], se.la[both.genes,i], method="pearson")
  plot(sm.la[both.genes,i], se.la[both.genes,i], main=paste0(i, " (",round(c1,2),")"),
       xlab="PRECAST logFC", ylab="Seurat logFC")
}




## L-R ----
# the layer-restricted model is asking fundamentally the different questions because the clusters are different
sm.lr = read.csv("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_smoothed-k9-1663_rev-gene-input_F-test_all.csv")
se.lr = read.csv("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_seurat-pc30_rev-gene-input_F-test_all.csv")
#cc.lr = read.csv("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_custom-cluster_rev-gene-input_F-test_all.csv")

both.genes = intersect(sm.lr$gene_id, se.lr$gene_id)
length(both.genes) #13015
rownames(sm.lr) <- sm.lr$gene_id
rownames(se.lr) <- se.lr$gene_id

par(mfrow=c(1,1))
plot(sm.lr[both.genes,"F"], se.lr[both.genes,"F"], xlab="PRECAST F statistic", ylab="Seurat F statistic",
     main="Layer-restricted", sub=paste("Spearman rho =",round(cor(sm.lr[both.genes,"F"], se.lr[both.genes,"F"], method="pearson"),3)))
cor(sm.lr[both.genes,"F"], se.lr[both.genes,"F"], method="pearson") # 0.6827596

# logFC
plist <- lapply(comparisons, function(i) {
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
    mutate(clusters=factor(clusters, levels=comp_names[[2]]), 
           PRECAST=factor(PRECAST, levels=rev(comp_names[[1]])))
  ggplot(cor.df, aes(x=clusters, y=PRECAST, fill=pearson_r))+
    geom_tile(color="grey50", linewidth=.1)+
    scale_fill_gradientn("Pearson\nrho", 
                         colors=colorRampPalette(RColorBrewer::brewer.pal(n=7, "RdBu")[7:1])(20),
                         limits=c(0, 1))+
    #geom_text(data=filter(cor.df, pearson_r>.5), aes(label=round(pearson_r,2)))+
    labs(title=paste(i, "L-R logFC"))+
    theme_minimal()+theme(aspect.ratio=1,
                          legend.key.width=unit(8,"pt"), legend.key.height=unit(10,"pt"), #legend.position="bottom",
                          legend.title = element_text(size=7), legend.text = element_text(size=6),
                          text=element_text(size=10),
                          axis.title.x=element_blank(), axis.title.y=element_blank(),
                          axis.ticks = element_line(color="grey50", linewidth=.3), axis.text.y=element_text(color="black"),
                          axis.text.x= element_text(angle=90, hjust=1, vjust=.5, size=6, color="black"))
})

do.call(gridExtra::grid.arrange, c(plist, ncol=2))

#nested F test
sm.lr = read.csv("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_smoothed-k9-1663_rev-gene-input_F-test_each.csv") %>%
  mutate(cluster=factor(cluster, levels=comp_names[[1]]))
se.lr = read.csv("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_seurat-pc30_rev-gene-input_F-test_each.csv") %>%
  mutate(cluster=factor(cluster, levels=comp_names[[2]]))

sm.wide = tidyr::pivot_wider(sm.lr[,c("gene_id","cluster","F")], names_from="cluster", values_from="F")
se.wide = tidyr::pivot_wider(se.lr[,c("gene_id","cluster","F")], names_from="cluster", values_from="F")

both.wide = left_join(filter(sm.wide, gene_id %in% both.genes),
                      filter(se.wide, gene_id %in% both.genes), by=c("gene_id"),
                      suffix=c("_sm","_se"))
head(both.wide)
wide.mtx = as.matrix(both.wide[,-1])
rownames(wide.mtx) <- both.wide$gene_id

corr.m = cor(wide.mtx, method="pearson")
dim(corr.m)
corr.m2 = corr.m[1:6, 7:ncol(corr.m)]
rownames(corr.m2) <- comp_names[[1]]
colnames(corr.m2) <- comp_names[[2]]


cor.df = tidyr::pivot_longer(tibble::rownames_to_column(as.data.frame(corr.m2), var="PRECAST"), all_of(comp_names[[2]]),
                    names_to="clusters", values_to="pearson_r") %>%
  mutate(clusters=factor(clusters, levels=comp_names[[2]]), 
         PRECAST=factor(PRECAST, levels=rev(comp_names[[1]])))

ggplot(cor.df, aes(x=clusters, y=PRECAST, fill=pearson_r))+
  geom_tile(color="grey50", linewidth=.1)+
  scale_fill_gradientn("Pearson\nrho", 
                       colors=colorRampPalette(RColorBrewer::brewer.pal(n=7, "RdBu")[7:1])(20),
                       limits=c(-1, 1))+#limits=c(0, 1))+
  #geom_text(data=filter(cor.df, pearson_r>.5), aes(label=round(pearson_r,2)))+
  labs(title="F statistic (L-R nested)")+
  theme_minimal()+theme(aspect.ratio=1,
                        legend.key.width=unit(8,"pt"), legend.key.height=unit(10,"pt"), #legend.position="bottom",
                        legend.title = element_text(size=7), legend.text = element_text(size=6),
                        text=element_text(size=10),
                        axis.title.x=element_blank(), axis.title.y=element_blank(),
                        axis.ticks = element_line(color="grey50", linewidth=.3), axis.text.y=element_text(color="black"),
                        axis.text.x= element_text(angle=90, hjust=1, vjust=.5, color="black"))


