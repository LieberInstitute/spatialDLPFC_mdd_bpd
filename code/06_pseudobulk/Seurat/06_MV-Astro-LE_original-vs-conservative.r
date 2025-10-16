setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(dplyr)
	library(ggplot2)
	library(UpSetR)
})
set.seed(123)

cpList <- readRDS("plots/colorPalettes.rds")
source("code/06_pseudobulk/custom_functions.r")

orig_all = read.csv("processed-data/06_pseudobulk/Seurat/layer-enrichment_seurat-pc30_all-results.csv")
cons_all = read.csv("processed-data/06_pseudobulk/Seurat/layer-enrichment_conservative_seurat-pc30_all-results.csv")
szbd_all = read.csv("processed-data/06_pseudobulk/SZBDMulti-seq/layer-enrichment_control-low-res_all-results.csv")

upsetList = lapply(c("Micro.Vasc","Astro"), function(x) {
  orig_large = filter(orig_all, seurat_label==x, sig_group=="large effect")$gene_id
  cons_large = filter(cons_all, seurat_label==x, sig_group=="large effect")$gene_id
  szbd_large = filter(szbd_all, seurat_label==x, sig_group=="large effect")
  szbd_large.other = filter(szbd_all, seurat_label!=x, sig_group=="large effect")
  szbd_large = setdiff(szbd_large$gene_id, szbd_large.other$gene_id)
  
  largeList = list("Orig"=orig_large, 
                   "Cons"=cons_large,
                   "SZBD"=szbd_large)
  return(largeList)
})
names(upsetList) = c("Micro.Vasc","Astro")

ilist = list(list("Cons"), list("Orig"), list("SZBD"),
             list("SZBD","Orig"), list("SZBD","Cons"), list("Orig","Cons"),
             list("SZBD","Orig","Cons"))


#dotplots
#load("processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo-dotplot_azimuth-broad.Rdata")
#colnames(sce_summ) <- gsub(" ", "_", gsub("/", "\\.", sce_summ$azimuth_broad))
load("processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo-dotplot_seurat-low-res.Rdata")
colnames(sce_summ) <- sce_summ$seurat_low.res
plot.genes = c(filter(szbd_all, seurat_label=="Micro.Vasc", adj.P.Val<1e-10) %>% 
                  slice_max(n=10, logFC) %>% pull(gene_id),
                filter(szbd_all, seurat_label=="Astro", adj.P.Val<1e-10) %>% 
                  slice_max(n=10, logFC) %>% pull(gene_id))
names(plot.genes) = rowData(sce_summ)[plot.genes,"gene_name"]

dot.df = dotplotDF(sce_summ, plot.genes, row_data=NULL)
dot.df$gene_name = rowData(sce_summ)[dot.df$gene_id,"gene_name"]
dot.df$gene_name_f = factor(dot.df$gene_name, levels=rev(names(plot.genes)))

dot.df$clusters = factor(dot.df$clusters, levels=levels(sce_summ$seurat_low.res))
dot.df$fill_color = factor(dot.df$gene_name, levels=names(plot.genes),
                            labels=c(rep("Micro.Vasc", 10), rep("Astro", 10)))

p1 <- ggplot(dot.df, aes(x=clusters, y=gene_name_f, 
                    color=mean_expr_scaled, size=prop_spots))+
  geom_tile(aes(fill=fill_color, size=NA), alpha=.5, color="transparent", show.legend=c(size=FALSE))+
  scale_fill_manual("Top 10 SZBD\nLE genes", values=cpList$transfer.light)+
  geom_count()+scale_color_gradient("Avg. expr.\n(scaled)", low="white", high="black")+
  scale_size("Prop. of\nspots", range=c(1,6), limits=c(0,1), breaks=c(0,.5,1))+
  labs(title="Top 10 M.V and Astro from SZBD control (low.res)", x="SZBD control pseudobulk")+
  theme_minimal()+theme(axis.title.y=element_blank(), axis.text.y=element_text(face="italic"))


#original
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-dotplot_dx-sex-seurat-pc30.Rdata")
cond_sex = c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M")
seurat_levels= c("Micro.Vasc","Astro","L2.3","L4","L5","L6","Oligo","Inhb")
spe_summ$seurat_label = factor(spe_summ$seurat_label, levels=seurat_levels)
spe_summ$sample_id = factor(paste(spe_summ$condition, spe_summ$sex, spe_summ$seurat_label),
                            levels=as.character(outer(cond_sex, seurat_levels, paste)))
colnames(spe_summ) <- spe_summ$sample_id
spe_summ_orig <- spe_summ

dot.df_orig = dotplotDF(spe_summ_orig, plot.genes, summarize_groups=T, cluster_labels="seurat_label", row_data=NULL)

dot.df_orig$gene_name = rowData(spe_summ_orig)[dot.df_orig$gene_id,"gene_name"]
dot.df_orig$gene_name_f = factor(dot.df_orig$gene_name, levels=rev(names(plot.genes)))

dot.df_orig$clusters = factor(dot.df_orig$clusters, levels=seurat_levels)
dot.df_orig$fill_color = factor(dot.df_orig$gene_name, levels=names(plot.genes),
                                 labels=c(rep("Micro.Vasc", 10), rep("Astro", 10)))

p2 <- ggplot(dot.df_orig, aes(x=clusters, y=gene_name_f, 
                         color=mean_expr_scaled, size=prop_spots))+
  geom_tile(aes(fill=fill_color, size=NA), alpha=.5, color="transparent", show.legend=c(size=FALSE))+
  scale_fill_manual("Top 10 SZBD\nLE genes", values=cpList$transfer.light)+
  geom_count()+scale_color_gradient("Avg. expr.\n(scaled)", low="white", high="black")+
  scale_size("Prop. of\nspots", range=c(1,6), limits=c(0,1), breaks=c(0,.5,1))+
  labs(title="Top 10 M.V and Astro from SZBD control (low.res)", x="Seurat PC30 (original) pseudobulk")+
  theme_minimal()+theme(axis.title.y=element_blank(), axis.text.y=element_text(face="italic"))


#conservative
load("processed-data/06_pseudobulk/Seurat/spe_n119_conservative_pseudo-dotplot_dx-sex-seurat-pc30.Rdata")
cond_sex = c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M")
seurat_levels= c("Micro.Vasc","Astro","L2.3","L4","L5","L6","Oligo","Inhb")
spe_summ$seurat_label = factor(spe_summ$seurat_label, levels=seurat_levels)
spe_summ$sample_id = factor(paste(spe_summ$condition, spe_summ$sex, spe_summ$seurat_label),
                            levels=as.character(outer(cond_sex, seurat_levels, paste)))
colnames(spe_summ) <- spe_summ$sample_id
spe_summ_cons <- spe_summ

dot.df_cons = dotplotDF(spe_summ_cons, plot.genes, summarize_groups=T, cluster_labels="seurat_label", row_data=NULL)

dot.df_cons$gene_name = rowData(spe_summ_cons)[dot.df_cons$gene_id,"gene_name"]
dot.df_cons$gene_name_f = factor(dot.df_cons$gene_name, levels=rev(names(plot.genes)))

dot.df_cons$clusters = factor(dot.df_cons$clusters, levels=seurat_levels)
dot.df_cons$fill_color = factor(dot.df_cons$gene_name, levels=names(plot.genes),
                                labels=c(rep("Micro.Vasc", 10), rep("Astro", 10)))

p3 <- ggplot(dot.df_cons, aes(x=clusters, y=gene_name_f, 
                              color=mean_expr_scaled, size=prop_spots))+
  geom_tile(aes(fill=fill_color, size=NA), alpha=.5, color="transparent", show.legend=c(size=FALSE))+
  scale_fill_manual("Top 10 SZBD\nLE genes", values=cpList$transfer.light)+
  geom_count()+scale_color_gradient("Avg. expr.\n(scaled)", low="white", high="black")+
  scale_size("Prop. of\nspots", range=c(1,6), limits=c(0,1), breaks=c(0,.5,1))+
  labs(title="Top 10 M.V and Astro from SZBD control (low.res)", x="Seurat PC30 (conservative) pseudobulk")+
  theme_minimal()+theme(axis.title.y=element_blank(), axis.text.y=element_text(face="italic"))

pdf(file="plots/06_pseudobulk/Seurat/compare-MV-Astro_pc30-original-vs-pc30-conservative.pdf")
upset(fromList(upsetList[[1]]), intersections= ilist, keep.order = T, text.scale=2, mb.ratio=c(.5,.5), 
      sets.x.label = "Micro.Vasc\nSet Size", mainbar.y.label = "Micro.Vasc\nIntersection Size")
upset(fromList(upsetList[[2]]), intersections= ilist, keep.order = T, text.scale=2, mb.ratio=c(.5,.5),
      sets.x.label = "Astro Set Size", mainbar.y.label = "Astro\nIntersection Size")
p1
p2
p3
dev.off()
cat("\n\nPDF saved to: plots/06_pseudobulk/Seurat/compare-MV-Astro_pc30-original-vs-pc30-conservative.pdf\n")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
