setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(edgeR)
	library(dplyr)
	library(ggplot2)
	library(scater)
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
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-no-lowUMI-heatmap_dx-sex-seurat-pc30.Rdata")
cond_sex = c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M")
seurat_levels= c("Micro.Vasc","Astro","L2.3","L4","L5","L6","Oligo","Inhb")
spe_summ$sample_id = factor(paste(spe_summ$condition, spe_summ$sex, spe_summ$seurat_qual.genes_pc30.kweight50),
                                 levels=as.character(outer(cond_sex, seurat_levels, paste)))
colnames(spe_summ) <- spe_summ$sample_id

#load in enrichment results and format
resList <- readRDS("processed-data/06_pseudobulk/Seurat/lmFit-list_seurat-pc30-no-lowUMI_covars-condition-sex-nspots-pc3.rda")
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

t_pc30 <- do.call(cbind, lapply(enrichList, function(x) x[,c(grep("^t_", colnames(x), value=T),"gene_id","gene_name")]))
t_pc30 <- t_pc30[,grep("^t_", colnames(t_pc30))]
colnames(t_pc30) <- gsub("^t_", "", colnames(t_pc30))

lf_pc30 <- do.call(cbind, lapply(enrichList, function(x) x[,c(grep("^logFC_", colnames(x), value=T),"gene_id","gene_name")]))
lf_pc30 <- lf_pc30[,grep("^logFC_", colnames(lf_pc30))]
colnames(lf_pc30) <- gsub("^logFC_", "", colnames(lf_pc30))


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
enrich.df$sig_group = ifelse(enrich.df$logFC>1, "large effect", enrich.df$adj.P.Val_bin2)
table(enrich.df[,c("sig_group","adj.P.Val_bin2")], useNA="ifany")
cat("\n\n")

enrich.df$seurat_label = factor(enrich.df$seurat_label, levels=c("Micro.Vasc","Astro","L2.3","L4","L5","L6","Oligo","Inhb"))
enrich.df_pc30 = left_join(enrich.df, rdata[,c("gene_id","gene_type")], by=c("gene_id")) %>%
  mutate(gene_type_ptn = gene_type=="protein_coding")
spe_summ_pc30 <- spe_summ

#save results
write.csv(t_pc30, "processed-data/06_pseudobulk/Seurat/layer-enrichment_seurat-pc30-no-lowUMI_t-stat.csv", row.names=T)
write.csv(lf_pc30, "processed-data/06_pseudobulk/Seurat/layer-enrichment_seurat-pc30-no-lowUMI_logFC.csv", row.names=T)
write.csv(enrich.df_pc30, "processed-data/06_pseudobulk/Seurat/layer-enrichment_seurat-pc30-no-lowUMI_all-results.csv", row.names=F)

cat("\n\nLayer enrichment results saved to:",
	"\n>> processed-data/06_pseudobulk/Seurat/layer-enrichment_seurat-pc30-no-lowUMI_t-stat.csv",
	"\n>> processed-data/06_pseudobulk/Seurat/layer-enrichment_seurat-pc30-no-lowUMI_logFC.csv",
	"\n>> processed-data/06_pseudobulk/Seurat/layer-enrichment_seurat-pc30-no-lowUMI_all-results.csv\n\n")


#make plots
cpList <- readRDS("plots/colorPalettes.rds")

#load in dlpfc marker genes (manually curated)
source("code/06_pseudobulk/dlpfc_genes.r")

dlpfc.genes.ids = rownames(rdata)[rdata$gene_name %in% unlist(dlpfc.genes)]
cat("\nNumber of dlPFC genes to highlight:", length(dlpfc.genes.ids),"\n\n")

#load in dotplotDF function and format dataframe for dotplot
source("code/06_pseudobulk/custom_functions.r")
pc30.df = dotplotDF(spe_summ_pc30, dlpfc.genes.ids, summarize_groups=T, 
                    cluster_labels="seurat_qual.genes_pc30.kweight50") %>%
  mutate(gene_name_f=factor(gene_name, levels=rev(unlist(dlpfc.genes))))

#dlpfc marker dotplot
c1 = unlist(lapply(names(dlpfc.genes), function(x) rep(x, length.out=length(dlpfc.genes[[x]]))))
pc30.df2 = mutate(pc30.df, fill_color = factor(gene_name, levels=unlist(dlpfc.genes), labels=c1))
p1 <- ggplot(pc30.df, aes(x=clusters, y=gene_name_f, color=mean_expr_scaled, size=prop_spots))+
  geom_tile(data=pc30.df2, aes(fill=fill_color, size=NA), alpha=.5, color="transparent", show.legend=c(size=FALSE))+
  scale_fill_manual(values=cpList$low.res.light, guide="none")+
  geom_count()+scale_color_gradient(low="white", high="black")+
  scale_x_discrete(labels=c("Micro\nVasc","Astro","L2.3","L4","L5","L6","Oligo","Inhb"))+
  scale_size(range=c(1,6), limits=c(0,1), breaks=c(0,.5,1))+
  labs(color="Avg. expr.\n(scaled)", size="Prop. of\nspots",
       y="common marker genes", title="MBv label transfer")+
  theme_minimal()+theme(axis.title.x=element_blank(), 
                        axis.text.y=element_text(face="italic"), legend.key.size=unit(15,"pt"),
                        axis.title.y=element_text(margin=margin(0,20,0,20,"pt")))

#volcano plots
#by protein coding genes (with fixed y)
library(ggh4x)
strip <- strip_themed(background_x = elem_list_rect(fill=c("white","white")),
                      background_y = elem_list_rect(fill = cpList$transfer.light, alpha=.5))
p2 <- ggplot(enrich.df_pc30, aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point(size=.1)+
  facet_grid2(rows=vars(seurat_label), cols=vars(gene_type_ptn),
             labeller=as_labeller(c("FALSE"="Not protein coding","TRUE"="Protein coding",
                                    Micro.Vasc="Micro\nVasc",Astro="Astro",
                                    L2.3="L2.3",L4="L4",L5="L5",L6="L6",
                                    Oligo="Oligo",Inhb="Inhb")),
             strip=strip)+
  xlim(-5.6,5.6)+
  labs(title="MBv label transfer", subtitle="Fixed y-axis limits across all panels")+
  theme_bw()+theme(panel.grid.minor=element_blank(), 
                   strip.text.y=element_text(angle=0))


#label dlpfc markers in protein coding genes (free y)
pc30_map = dlpfc.genes
pc30_map[["L2"]] <- c(pc30_map[["L2"]], pc30_map[["L3"]])
pc30_map[["L3"]] <- NULL
names(pc30_map)[grep("L2", names(pc30_map))] <- "L2.3"

enrich.df_pc30$plot_genes = FALSE
for(i in names(pc30_map)) {
  enrich.df_pc30$plot_genes = ifelse((enrich.df_pc30$gene_name %in% pc30_map[[i]]) & 
                                       (enrich.df_pc30$seurat_label==i), T, enrich.df_pc30$plot_genes)
}

## this is for the other version of the labeled volcano plot (see p3 in SZBDMulti-seq plot_layer-enrichment code)
#enrich.df_pc30$gene_color = NA
#for(i in names(dlpfc.genes)) {
#  enrich.df_pc30$gene_color = ifelse((enrich.df_pc30$gene_name %in% dlpfc.genes[[i]]) & 
#                                       (enrich.df_pc30$plot_genes==T), i, enrich.df_pc30$gene_color)
#}

p4 <- ggplot(filter(enrich.df_pc30, gene_type_ptn==T), 
             aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point(data=filter(enrich.df_pc30, gene_type_ptn==T, plot_genes==F), color="grey70", size=.3)+
  geom_point(data=filter(enrich.df_pc30, gene_type_ptn==T, plot_genes==T), color="black", size=1)+
  facet_grid2(rows= vars(seurat_label), scales="free_y",
              labeller=as_labeller(c(Micro.Vasc="Micro\nVasc",Astro="Astro",
                                     L2.3="L2.3",L4="L4",L5="L5",L6="L6",
                                     Oligo="Oligo",Inhb="Inhb")),
              strip=strip)+
  ggrepel::geom_text_repel(data= filter(enrich.df_pc30, plot_genes==T),
                            aes(label=gene_name), 
                            size=3, fontface="bold.italic", 
                            min.segment.length=0, max.overlaps=Inf)+
  xlim(-5.6,5.6)+
  labs(title="MBv label transfer", subtitle="Protein coding genes only")+
  theme_bw()+theme(panel.grid.minor=element_blank(), strip.text.y=element_text(angle=0),
                   plot.margin = margin(.2,2.5,.2,2.5,"cm"))

pdf(file="plots/06_pseudobulk/Seurat/seurat-pc30-no-lowUMI_layer-enrichment_plots.pdf",
    width=6, height=8)
p2
p4
p1
dev.off()

cat("\nPlots saved to: plots/06_pseudobulk/Seurat/seurat-pc30-no-lowUMI_layer-enrichment_plots.pdf\n")


#top top layer markers
cat("\nTop layer markers (adj. p<1e-30, logFC>1\n")
filter(enrich.df_pc30, adj.P.Val<1e-30, sig_group=="large effect") %>%
  group_by(seurat_label) %>% tally()


load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-no-lowUMI_sample-seurat-pc30_norm-filt.Rdata")

top.list = lapply(levels(enrich.df_pc30$seurat_label), function(x) {
  tmp = filter(enrich.df_pc30, adj.P.Val<1e-30, sig_group=="large effect", seurat_label==x)$gene_name
  phm = plotGroupedHeatmap(spe_pseudo, features=tmp, swap_rownames="gene_name",
                     group="seurat_label", center=T, cluster_cols=F, angle_col=0, silent=T)
  phm[[4]]
  })
names(top.list) = levels(enrich.df_pc30$seurat_label)

phm.list = gridExtra::marrangeGrob(top.list, ncol=1, nrow=1, top=quote(names(top.list)[g]))
ggsave("plots/06_pseudobulk/Seurat/seurat-pc30-no-lowUMI_layer-markers-adjp-1e30-logfc-1_heatmap.pdf", phm.list, width=7, height=11)

cat("\n\nTop layer markers (adj p<1e-30, logFC>1) heatmap saved to: plots/06_pseudobulk/Seurat/seurat-pc30-no-lowUMI_layer-markers-adjp-1e30-logfc-1_heatmap.pdf\n")



cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
