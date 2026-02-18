setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
})
set.seed(123)

source("code/06_pseudobulk/custom_functions.r")
cpList <- readRDS("plots/colorPalettes.rds")
cpList$low.res.bright
c("Micro.Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb")

#load in enrichment results
enrich.df_sn <- read.csv("processed-data/06_pseudobulk/SZBDMulti-seq/layer-enrichment_control-low-res_all-results.csv")
enrich.df_sn$seurat_label = factor(enrich.df_sn$seurat_label, levels=names(cpList$low.res.bright))
t_sn <- read.csv("processed-data/06_pseudobulk/SZBDMulti-seq/layer-enrichment_control-low-res_t-stat.csv", row.names=1)
lf_sn <- read.csv("processed-data/06_pseudobulk/SZBDMulti-seq/layer-enrichment_control-low-res_logFC.csv", row.names=1)

enrich.df_smooth <- read.csv("processed-data/06_pseudobulk/PRECAST_smoothed/layer-enrichment_smoothed-k9-1663_all-results.csv")
enrich.df_smooth$smoothed = factor(enrich.df_smooth$smoothed, levels=names(cpList$smoothed.bright))
t_sm <-	read.csv("processed-data/06_pseudobulk/PRECAST_smoothed/layer-enrichment_smoothed-k9-1663_t-stat.csv", row.names=1)
lf_sm <- read.csv("processed-data/06_pseudobulk/PRECAST_smoothed/layer-enrichment_smoothed-k9-1663_logFC.csv", row.names=1)

#### https://www.bioconductor.org/packages/release/data/experiment/vignettes/spatialLIBD/inst/doc/spatialLIBD.html#spatiallibd-functions
#### We already covered fetch_data() which allows you to download the Human DLPFC Visium data from LIBD researchers and colleagues (Maynard, Collado-Torres, Weber et al., 2021).
### in the user guide/ function list I found the key for the 2024 dlPFC results
layer_modeling_results <- spatialLIBD::fetch_data(type = "spatialDLPFC_Visium_modeling_results")
t1 = layer_modeling_results$enrichment[,c(grep("t_stat",colnames(layer_modeling_results$enrichment), value=T),"ensembl")]
#colnames(t1) = gsub("t_stat_","", colnames(t1))
#colnames(t1) = gsub("Layer","L", colnames(t1))
#m1 = as.matrix(t1[,1:7])
fix.names = c("Sp09D01"="Mng","Sp09D02"="L1","Sp09D03"="L2","Sp09D04"="L5","Sp09D05"="L3","Sp09D06"="WM.1","Sp09D07"="L6","Sp09D08"="L4","Sp09D09"="WM.2")
m1 = as.matrix(t1[,1:9])
colnames(m1) <- fix.names
rownames(m1) = t1$ensembl
dim(m1) #12225 9


top.t = getTopGenes(t_sm, top_n=50)
top.lf = getTopGenes(lf_sm, top_n=50)
top50.both= union(top.t, top.lf)


#top 50 t stat
top100.list <- list("SZBDMulti-seq"=getTopGenes(t_sn, top_n=100),
                   "PRECAST (smoothed)"=getTopGenes(t_sm, top_n=100),
                   "spatialDLPFC"=getTopGenes(m1, top_n=100)
)
sapply(top100.list, length)

orderList <- list("spatialDLPFC"=c("Mng","L1","L2","L3","L4","L5","L6","WM.1","WM.2"),
                  "SZBDMulti-seq"=c("Micro.Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"),
                  "PRECAST (smoothed)"=c("L1","L2","L3.4","L5","L6","WM"))




both.genes1 = union(intersect(top100.list[["PRECAST (smoothed)"]], rownames(t_sn)),
                        intersect(top100.list[["SZBDMulti-seq"]], rownames(t_sm)))
length(both.genes1) #1253
table(rdata[both.genes1,"gene_type"])
#lncRNA protein_coding 
#390            863
table(rdata[intersect(both.genes1,top100.list[["PRECAST (smoothed)"]]),"gene_type"])
#lncRNA protein_coding 
#67            514
round(67/(514+67),3)
table(rdata[intersect(both.genes1,top100.list[["SZBDMulti-seq"]]),"gene_type"])
#lncRNA protein_coding 
#362            452 
round(362/(452+362),3)


cor.res = cor(t_sn[both.genes1,],t_sm[both.genes1,])
cor.df1 = tibble::rownames_to_column(as.data.frame(cor.res), var="clusters") %>%
  tidyr::pivot_longer(all_of(colnames(cor.res)), names_to="PRECAST", 
                      values_to="pearson_r")
cor.df1 = mutate(cor.df1, clusters=factor(clusters, levels=orderList[["SZBDMulti-seq"]], labels=c("M.V","Astro","L2","L3","L4","L5","L6","Oligo","Inhb")),
       compare_to="SZBDMulti-seq")


both.genes2 = union(intersect(top100.list[["PRECAST (smoothed)"]], rownames(m1)),
                    intersect(top100.list[["spatialDLPFC"]], rownames(t_sm)))
length(both.genes2) #942

cor.res = cor(m1[both.genes2,],t_sm[both.genes2,])
cor.df2 = tibble::rownames_to_column(as.data.frame(cor.res), var="clusters") %>%
  tidyr::pivot_longer(all_of(colnames(cor.res)), names_to="PRECAST", 
                      values_to="pearson_r")
cor.df2 = mutate(cor.df2, clusters=factor(clusters, levels=orderList[["spatialDLPFC"]]),
                 compare_to="spatialDLPFC")
head(cor.df1)
head(cor.df2)

cor.df = bind_rows(cor.df1, cor.df2) %>%
  mutate(clusters= factor(clusters, levels=c("Mng","M.V","L1","Astro","L2","L3","L4","L5","L6","Oligo","WM.1","WM.2","Inhb")),
         PRECAST= factor(PRECAST, levels=rev(orderList[["PRECAST (smoothed)"]])),
         compare_to= factor(compare_to, levels=c("spatialDLPFC","SZBDMulti-seq")))

p1 <- ggplot(cor.df, aes(x=clusters, y=PRECAST, fill=pearson_r))+
  geom_tile(color="grey50", linewidth=.1)+
  scale_fill_gradientn("Pearson\nrho", 
                       colors=colorRampPalette(RColorBrewer::brewer.pal(n=7, "RdBu")[7:1])(20),
                       limits=c(-1, 1))+
  facet_grid(cols=vars(compare_to), scales="free_x")+
  theme_minimal()+theme(panel.grid= element_blank(), aspect.ratio=1,
                        legend.key.width=unit(8,"pt"), legend.key.height=unit(10,"pt"), #legend.position="bottom",
                        legend.title = element_text(size=7), legend.text = element_text(size=6),
                        text=element_text(size=10),
                        axis.title.x=element_blank(), axis.title.y=element_blank(),
                        axis.ticks = element_line(color="grey50", linewidth=.3), axis.text.y=element_text(color="black"),
                        axis.text.x= element_text(angle=90, hjust=1, vjust=.5, size=6, color="black"))


ggsave(file="plots/publication/Figure1/spatial-registration.pdf",
       p1,
       width=3, height=2)


length(both.genes1) #1253
r1 = rdata[both.genes1,]
table(r1$gene_type)
#lncRNA protein_coding 
#390            863 
r2 = r1[r1$gene_type=="protein_coding","gene_id"]
length(r2) #863

cor.res = cor(t_sn[r2,],t_sm[r2,])
cor.df3 = tibble::rownames_to_column(as.data.frame(cor.res), var="clusters") %>%
  tidyr::pivot_longer(all_of(colnames(cor.res)), names_to="PRECAST", 
                      values_to="pearson_r")
cor.df3 = mutate(cor.df3, clusters=factor(clusters, levels=orderList[["SZBDMulti-seq"]], labels=c("M.V","Astro","L2","L3","L4","L5","L6","Oligo","Inhb")),
                 compare_to="SZBDMulti-seq")

length(both.genes2) #942
r3 = rdata[both.genes2,]
table(r3$gene_type)
#lncRNA protein_coding 
#32            910 
r4 = r3[r3$gene_type=="protein_coding","gene_id"]
length(r4) #910

cor.res = cor(m1[r4,],t_sm[r4,])
cor.df4 = tibble::rownames_to_column(as.data.frame(cor.res), var="clusters") %>%
  tidyr::pivot_longer(all_of(colnames(cor.res)), names_to="PRECAST", 
                      values_to="pearson_r")
cor.df4 = mutate(cor.df4, clusters=factor(clusters, levels=orderList[["spatialDLPFC"]]),
                 compare_to="spatialDLPFC")


cor.df = bind_rows(cor.df3, cor.df4) %>%
  mutate(clusters= factor(clusters, levels=c("Mng","M.V","L1","Astro","L2","L3","L4","L5","L6","Oligo","WM.1","WM.2","Inhb")),
         PRECAST= factor(PRECAST, levels=rev(orderList[["PRECAST (smoothed)"]])),
         compare_to= factor(compare_to, levels=c("spatialDLPFC","SZBDMulti-seq")))

p2 <- ggplot(cor.df, aes(x=clusters, y=PRECAST, fill=pearson_r))+
  geom_tile(color="grey50", linewidth=.1)+
  scale_fill_gradientn("Pearson\nrho", 
                       colors=colorRampPalette(RColorBrewer::brewer.pal(n=7, "RdBu")[7:1])(20),
                       limits=c(-1, 1))+
  facet_grid(cols=vars(compare_to), scales="free_x")+
  theme_minimal()+theme(panel.grid= element_blank(), aspect.ratio=1,
                        legend.key.width=unit(8,"pt"), legend.key.height=unit(10,"pt"), #legend.position="bottom",
                        legend.title = element_text(size=7), legend.text = element_text(size=6),
                        text=element_text(size=10),
                        axis.title.x=element_blank(), axis.title.y=element_blank(),
                        axis.ticks = element_line(color="grey50", linewidth=.3), axis.text.y=element_text(color="black"),
                        axis.text.x= element_text(angle=90, hjust=1, vjust=.5, size=6, color="black"))

ggsave(file="plots/publication/Figure1/spatial-registration_protein-coding-only.pdf",
       p2,
       width=3, height=2)
