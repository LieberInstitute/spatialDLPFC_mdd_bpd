setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(edgeR)
	library(dplyr)
	library(ggplot2)
})
set.seed(123)

# load in unfiltered spe object which should have whole gene universe
### the information here (https://www.synapse.org/Synapse:syn22963646) indicates that they used
### the hg38 genome for their reference which is what we used too (?)
spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_")
rdata = as.data.frame(rowData(spe))
dim(rdata) #36601 7
rm(spe)

#load in sce for heatmap/dotplots
load("processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo-heatmap_seurat-low-res.Rdata")
colnames(sce_summ) <- sce_summ$seurat_low.res

#load in enrichment results and format
resList <- readRDS("processed-data/06_pseudobulk/SZBDMulti-seq/lmFit-list_control-low-res_covars-sex-ncells.rda")
enrichList = lapply(names(resList), function(x) {
  tmp = resList[[x]]
  tmp = eBayes(tmp)
  tmp = topTable(tmp, "res", #p.value=.05, 
                 n=Inf, sort.by="none")
  colnames(tmp) = paste(colnames(tmp), x, sep="_")
  tmp$gene_id = rownames(tmp)
  tmp$gene_name = rdata[rownames(tmp),"gene_name"]
  return(tmp)
})

t_sn <- do.call(cbind, lapply(enrichList, function(x) x[,c(grep("^t_", colnames(x), value=T),"gene_id","gene_name")]))
t_sn <- t_sn[,grep("^t_", colnames(t_sn))]
colnames(t_sn) <- gsub("^t_", "", colnames(t_sn))

lf_sn <- do.call(cbind, lapply(enrichList, function(x) x[,c(grep("^logFC_", colnames(x), value=T),"gene_id","gene_name")]))
lf_sn <- lf_sn[,grep("^logFC_", colnames(lf_sn))]
colnames(lf_sn) <- gsub("^logFC_", "", colnames(lf_sn))

#compile results in long form dframe
enrich.df = do.call(rbind, lapply(enrichList, function(x) {
  c1 = strsplit(colnames(x)[1:6], "_")
  cell.type = c1[[1]][2]
  colnames(x)[1:6] = sapply(c1, function(y) y[[1]])
  x$seurat_label = cell.type
  rownames(x) <- NULL
  return(x)
}))
enrich.df$adj.P.Val_bin = cut(enrich.df$adj.P.Val, breaks=c(0,1e-10,.05,1), include.lowest=T)
table(enrich.df$adj.P.Val_bin, useNA="ifany")
cat("\n\n")
enrich.df$adj.P.Val_bin2 = as.character(factor(enrich.df$adj.P.Val_bin, levels=levels(enrich.df$adj.P.Val_bin), labels=c("highly sig.","sig.","NS")))
enrich.df$sig_group = ifelse(enrich.df$logFC>2, "large effect", enrich.df$adj.P.Val_bin2)
table(enrich.df[,c("sig_group","adj.P.Val_bin2")], useNA="ifany")
cat("\n\n")
enrich.df$sig_group2 = ifelse(enrich.df$sig_group=="large effect" & enrich.df$adj.P.Val_bin2=="highly sig.", "large effect", enrich.df$adj.P.Val_bin2)
table(enrich.df[,c("sig_group2","adj.P.Val_bin2")], useNA="ifany")
cat("\n\n")

enrich.df$seurat_label = factor(enrich.df$seurat_label, levels=c("Micro.Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"))
enrich.df_sn <- left_join(enrich.df, rdata[,c("gene_id","gene_type")], by=c("gene_id")) %>%
  mutate(gene_type_ptn = gene_type=="protein_coding")

#save results
write.csv(t_sn, "processed-data/06_pseudobulk/SZBDMulti-seq/layer-enrichment_control-low-res_t-stat.csv", row.names=T)
write.csv(lf_sn, "processed-data/06_pseudobulk/SZBDMulti-seq/layer-enrichment_control-low-res_logFC.csv", row.names=T)
write.csv(enrich.df_sn, "processed-data/06_pseudobulk/SZBDMulti-seq/layer-enrichment_control-low-res_all-results.csv", row.names=F)

cat("\n\nLayer enrichment results saved to:",
	"\n>> processed-data/06_pseudobulk/SZBMulti-seq/layer-enrichment_control-low-res_t-stat.csv",
	"\n>> processed-data/06_pseudobulk/SZBMulti-seq/layer-enrichment_control-low-res_logFC.csv",
	"\n>> processed-data/06_pseudobulk/SZBMulti-seq/layer-enrichment_control-low-res_all-results.csv\n\n")


#make plots
cpList <- readRDS("plots/colorPalettes.rds")

#load in dlpfc marker genes (manually curated)
source("code/06_pseudobulk/dlpfc_genes.r")

dlpfc.genes.ids = rownames(rdata)[rdata$gene_name %in% unlist(dlpfc.genes)]
cat("\nNumber of dlPFC genes to highlight:", length(dlpfc.genes.ids),"\n\n")

#load in dotplotDF function and format dataframe for dotplot
source("code/06_pseudobulk/custom_functions.r")
sn.df = dotplotDF(sce_summ, dlpfc.genes.ids) %>% 
  mutate(clusters= factor(clusters, levels=levels(sce_summ$seurat_low.res)),
         gene_name_f=factor(gene_name, levels=rev(unlist(dlpfc.genes))))

#dlpfc marker dotplot
c1 = unlist(lapply(names(dlpfc.genes), function(x) rep(x, length.out=length(dlpfc.genes[[x]]))))
sn.df2 = mutate(sn.df, fill_color = factor(gene_name, levels=unlist(dlpfc.genes), labels=c1))
p1 <- ggplot(sn.df, aes(x=clusters, y=gene_name_f, color=mean_expr_scaled, size=prop_spots))+
  geom_tile(data=sn.df2, aes(fill=fill_color), alpha=.5, color="transparent", show.legend=c(size=FALSE))+
  scale_fill_manual("Gene\nmarker\nfor:", values=cpList$low.res.light, guide="none")+
  geom_count()+scale_color_gradient(low="white", high="black")+
  scale_x_discrete(labels=c("Micro\nVasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"))+
  scale_size(range=c(1,6), limits=c(0,1), breaks=c(0,.5,1))+
  labs(color="Avg. expr.\n(scaled)", size="Prop. of\nnuclei",
       y="common marker genes", title="SZBDMulti-seq")+
  theme_minimal()+theme(axis.title.x=element_blank(), 
                        axis.text.y=element_text(face="italic"), legend.key.size=unit(15,"pt"),
                        axis.title.y=element_text(margin=margin(0,20,0,0,"pt")))


#volcano plots
#by protein coding genes (with fixed y)
library(ggh4x)
strip <- strip_themed(background_x = elem_list_rect(fill=c("white","white")),
                      background_y = elem_list_rect(fill = cpList$low.res.light, alpha=.5))
p2 <- ggplot(enrich.df_sn, aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point(size=.1)+
  facet_grid2(rows=vars(seurat_label), cols=vars(gene_type_ptn),
             labeller=as_labeller(c("FALSE"="Not protein coding","TRUE"="Protein coding",
                                    Micro.Vasc="Micro\nVasc",Astro="Astro",
                                    L2="L2",L3="L3",L4="L4",L5="L5",L6="L6",
                                    Oligo="Oligo",Inhb="Inhb")),
             strip=strip)+
  xlim(-10,10)+
  labs(title="SZBDMulti-seq", subtitle="Faceted by if gene is protein coding")+
  theme_bw()+theme(panel.grid.minor=element_blank(), 
                   #strip.background = element_rect(fill="transparent"),
                   strip.text.y=element_text(angle=0))

#label dlpfc markers in protein coding genes (free y)
enrich.df_sn$plot_genes = FALSE
for(i in names(dlpfc.genes)) {
  enrich.df_sn$plot_genes = ifelse((enrich.df_sn$gene_name %in% dlpfc.genes[[i]]) & 
                                     (enrich.df_sn$seurat_label==i), T, enrich.df_sn$plot_genes)
}
enrich.df_sn$gene_color = NA
for(i in names(dlpfc.genes)) {
  enrich.df_sn$gene_color = ifelse((enrich.df_sn$gene_name %in% dlpfc.genes[[i]]) & 
                                     (enrich.df_sn$plot_genes==T), i, enrich.df_sn$gene_color)
}


## version with colored label (decided I like the other version better but keeping code just in case)
#p3 <- ggplot(filter(enrich.df_sn, gene_type_ptn==T), 
#       aes(x=logFC, y=-log10(adj.P.Val)))+
#  geom_point(data=filter(enrich.df_sn, gene_type_ptn==T, plot_genes==F), color="grey50", size=.3)+
#  geom_point(data=filter(enrich.df_sn, gene_type_ptn==T, plot_genes==T), color="black", size=1)+
#  facet_grid(rows= vars(seurat_label), scales="free_y")+
#  ggrepel::geom_label_repel(data= filter(enrich.df_sn, plot_genes==T),
#                            aes(label=gene_name, fill=gene_color), 
#                            size=3, fontface="italic", 
#                            min.segment.length=0, max.overlaps=Inf)+
#  scale_fill_manual(values=cpList$low.res.light, guide="none")+
#  xlim(-10,10)+
#  labs(title="SZBDMulti-seq", subtitle="Protein coding genes only")+
#  theme_bw()+theme(panel.grid.minor=element_blank(), strip.text.y=element_text(angle=0),
#                   strip.background = element_rect(fill="transparent"),
#                   plot.margin = margin(.2,2,.2,2,"cm"))


p4 <- ggplot(filter(enrich.df_sn, gene_type_ptn==T), 
             aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point(data=filter(enrich.df_sn, gene_type_ptn==T, plot_genes==F), color="grey70", size=.3)+
  geom_point(data=filter(enrich.df_sn, gene_type_ptn==T, plot_genes==T), color="black", size=1)+
  facet_grid2(rows= vars(seurat_label), scales="free_y",
              labeller=as_labeller(c(Micro.Vasc="Micro\nVasc",Astro="Astro",
                                     L2="L2",L3="L3",L4="L4",L5="L5",L6="L6",
                                     Oligo="Oligo",Inhb="Inhb")),
              strip=strip)+
  ggrepel::geom_text_repel(data= filter(enrich.df_sn, plot_genes==T),
                            aes(label=gene_name), 
                            size=3, fontface="bold.italic", 
                            min.segment.length=0, max.overlaps=Inf)+
  xlim(-10,10)+
  labs(title="SZBDMulti-seq", subtitle="Protein coding genes only")+
  theme_bw()+theme(panel.grid.minor=element_blank(), strip.text.y=element_text(angle=0),
                   plot.margin = margin(.2,2.5,.2,2.5,"cm"))

pdf(file="plots/06_pseudobulk/SZBDMulti-seq/control-low-res_layer-enrichment_plots.pdf",
    width=6, height=8)
p2
p4
p1
dev.off()

cat("\nPlots saved to: plots/06_pseudobulk/SZBDMulti-seq/control-low-res_layer-enrichment_plots.pdf\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
