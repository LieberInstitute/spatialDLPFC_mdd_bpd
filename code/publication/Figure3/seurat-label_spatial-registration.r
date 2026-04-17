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

#load in enrichment results
enrich.df_sn <- read.csv("processed-data/06_pseudobulk/SZBDMulti-seq/layer-enrichment_control-low-res_all-results.csv")
enrich.df_sn$seurat_label = factor(enrich.df_sn$seurat_label, levels=names(cpList$low.res.bright))
t_sn <- read.csv("processed-data/06_pseudobulk/SZBDMulti-seq/layer-enrichment_control-low-res_t-stat.csv", row.names=1)
lf_sn <- read.csv("processed-data/06_pseudobulk/SZBDMulti-seq/layer-enrichment_control-low-res_logFC.csv", row.names=1)

enrich.df_smooth <- read.csv("processed-data/06_pseudobulk/Seurat/layer-enrichment_seurat-pc30_all-results.csv")
enrich.df_smooth$seurat_label = factor(enrich.df_smooth$seurat_label, levels=names(cpList$transfer.bright)[c(1:4,8,5:7)])
t_sm <-	read.csv("processed-data/06_pseudobulk/Seurat/layer-enrichment_seurat-pc30_t-stat.csv", row.names=1)
lf_sm <- read.csv("processed-data/06_pseudobulk/Seurat/layer-enrichment_seurat-pc30_logFC.csv", row.names=1)

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


top.t = getTopGenes(t_sm, top_n=100)
length(top.t)

orderList <- list("spatialDLPFC"=c("Mng","L1","L2","L3","L4","L5","L6","WM.1","WM.2"),
                  "SZBDMulti-seq"=c("Micro.Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"),
                  "Seurat label"=c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))


both.genes0 = intersect(top.t, rownames(t_sn))
length(both.genes0)

cor.res = cor(t_sn[both.genes0,],t_sm[both.genes0,])
cor.df0 = tibble::rownames_to_column(as.data.frame(cor.res), var="clusters") %>%
  tidyr::pivot_longer(all_of(colnames(cor.res)), names_to="seurat", 
                      values_to="pearson_r")
cor.df0 = mutate(cor.df0, clusters=factor(clusters, levels=orderList[["SZBDMulti-seq"]], labels=c("M.V","Astro","L2","L3","L4","Inhb","L5","L6","Oligo")),
                 compare_to="SZBDMulti-seq")


both.genes1 = intersect(top.t, rownames(m1))
length(both.genes1) #481

cor.res = cor(m1[both.genes1,],t_sm[both.genes1,])
cor.df1 = tibble::rownames_to_column(as.data.frame(cor.res), var="clusters") %>%
  tidyr::pivot_longer(all_of(colnames(cor.res)), names_to="seurat", 
                      values_to="pearson_r")
cor.df1 = mutate(cor.df1, clusters=factor(clusters, levels=orderList[["spatialDLPFC"]]),
                 compare_to="spatialDLPFC")


cor.df = bind_rows(cor.df0, cor.df1) %>%
  mutate(clusters= factor(clusters, levels=c("Mng","M.V","L1","Astro","L2","L3","L4","Inhb","L5","L6","Oligo","WM.1","WM.2")),
         seurat= factor(seurat, levels=rev(orderList[["Seurat label"]]), labels=rev(c("M.V","Ast","L2.3","L4","Inb","L5","L6","Olg"))),
         compare_to= factor(compare_to, levels=c("spatialDLPFC","SZBDMulti-seq")))

p0 <- ggplot(cor.df, aes(x=clusters, y=seurat, fill=pearson_r))+
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

ggsave(file="plots/publication/Figure3/seurat-label_spatial-registration_seurat-MBv-markers-only.pdf",
       p0,
       width=3, height=2)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
