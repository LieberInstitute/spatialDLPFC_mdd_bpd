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

#load in sce for heatmap and dotplots
load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo-dotplot_dx-sex-smoothed-n1663-k9.Rdata")
cond_sex = c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M")
precast_levels= c("L1","L2","L3.4","L5","L6","WM")
spe_summ$sample_id = factor(paste(spe_summ$condition, spe_summ$sex, spe_summ$smoothed_k9_1663),
                            levels=as.character(outer(cond_sex, precast_levels, paste)))
colnames(spe_summ) <- spe_summ$sample_id

#load in enrichment results and format them
resList <- readRDS("processed-data/06_pseudobulk/PRECAST_smoothed/lmFit-list_smoothed-k9-1663_covars-condition-sex-nspots-pc3.rda")
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

t_sm <- do.call(cbind, lapply(enrichList, function(x) x[,c(grep("^t_", colnames(x), value=T),"gene_id","gene_name")]))
t_sm <- t_sm[,grep("^t_", colnames(t_sm))]
colnames(t_sm) <- gsub("^t_", "", colnames(t_sm))

lf_sm <- do.call(cbind, lapply(enrichList, function(x) x[,c(grep("^logFC_", colnames(x), value=T),"gene_id","gene_name")]))
lf_sm <- lf_sm[,grep("^logFC_", colnames(lf_sm))]
colnames(lf_sm) <- gsub("^logFC_", "", colnames(lf_sm))

#compile results in long form dframe
enrich.df = do.call(rbind, lapply(enrichList, function(x) {
  c1 = strsplit(colnames(x)[1:6], "_")
  cell.type = c1[[1]][2]
  colnames(x)[1:6] = sapply(c1, function(y) y[[1]])
  x$smoothed = cell.type
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

enrich.df$smoothed = factor(enrich.df$smoothed, levels=c("L1","L2","L3.4","L5","L6","WM"))
enrich.df_smooth = left_join(enrich.df, rdata[,c("gene_id","gene_type")], by=c("gene_id")) %>%
  mutate(gene_type_ptn = gene_type=="protein_coding")
spe_summ_sm <- spe_summ

#save results
write.csv(t_sm, "processed-data/06_pseudobulk/PRECAST_smoothed/layer-enrichment_smoothed-k9-1663_t-stat.csv", row.names=T)
write.csv(lf_sm, "processed-data/06_pseudobulk/PRECAST_smoothed/layer-enrichment_smoothed-k9-1663_logFC.csv", row.names=T)
write.csv(enrich.df_smooth, "processed-data/06_pseudobulk/PRECAST_smoothed/layer-enrichment_smoothed-k9-1663_all-results.csv", row.names=F)

cat("\n\nLayer enrichment results saved to:",
	"\n>> processed-data/06_pseudobulk/PRECAST_smoothed/layer-enrichment_smoothed-k9-1663_t-stat.csv",
	"\n>> processed-data/06_pseudobulk/PRECAST_smoothed/layer-enrichment_smoothed-k9-1663_logFC.csv",
	"\n>> processed-data/06_pseudobulk/PRECAST_smoothed/layer-enrichment_smoothed-k9-1663_all-results.csv\n\n")


#make plots
cpList <- readRDS("plots/colorPalettes.rds")

#load in dlpfc marker genes (manually curated)
source("code/06_pseudobulk/dlpfc_genes.r")

dlpfc.genes.ids = rownames(rdata)[rdata$gene_name %in% unlist(dlpfc.genes)]
cat("\nNumber of dlPFC genes to highlight:", length(dlpfc.genes.ids),"\n\n")

#load in dotplotDF function and format dataframe for dotplot
source("code/06_pseudobulk/custom_functions.r")
sm.df = dotplotDF(spe_summ_sm, dlpfc.genes.ids, 
                  summarize_groups=T, cluster_labels="smoothed_k9_1663") %>%
  mutate(gene_name_f=factor(gene_name, levels=rev(unlist(dlpfc.genes))))

#dlpfc marker dotplot
c1 = unlist(lapply(names(dlpfc.genes), function(x) rep(x, length.out=length(dlpfc.genes[[x]]))))
sm.df2 = mutate(sm.df, fill_color = factor(gene_name, levels=unlist(dlpfc.genes), labels=c1))

p1 <- ggplot(sm.df, aes(x=clusters, y=gene_name_f, color=mean_expr_scaled, size=prop_spots))+
  geom_tile(data=sm.df2, aes(fill=fill_color), alpha=.5, color="transparent", show.legend=c(size=FALSE))+
  scale_fill_manual(values=cpList$low.res.light, guide="none")+
  geom_count()+scale_color_gradient(low="white", high="black")+
  scale_x_discrete(labels=c("L1\n","L2","L3.4","L5","L6","WM"))+
  scale_size(range=c(1,6), limits=c(0,1), breaks=c(0,.5,1))+
  labs(color="Avg. expr.\n(scaled)", size="Prop. of\nspots",
       y="common marker genes", title="PRECAST (smoothed)")+
  theme_minimal()+theme(axis.title.x=element_blank(), 
                        axis.text.y=element_text(face="italic"), legend.key.size=unit(15,"pt"),
                        axis.title.y=element_text(margin=margin(0,20,0,40,"pt")))


#volcano plots
#by protein coding genes (with fixed y)
library(ggh4x)
strip <- strip_themed(background_x = elem_list_rect(fill=c("white","white")),
                      background_y = elem_list_rect(fill = cpList$smoothed.light, alpha=.5))

p2 <- ggplot(enrich.df_smooth, aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point(size=.1)+
  facet_grid2(rows=vars(smoothed), cols=vars(gene_type_ptn),
             labeller=as_labeller(c("FALSE"="Not protein coding","TRUE"="Protein coding",
                                    L1="L1",
                                    L2="L2",L3.4="L3.4",L5="L5",L6="L6",
                                    WM="WM")),
             strip=strip)+
  xlim(-5,5)+
  labs(title="PRECAST (smoothed)", subtitle="Fixed y-axis limits across all panels")+
  theme_bw()+theme(panel.grid.minor=element_blank(), 
                   strip.text.y=element_text(angle=0))

#label dlpfc markers in protein coding genes (free y)
sm_map <- list(L1= c(dlpfc.genes[["Micro.Vasc"]], dlpfc.genes[["Astro"]]),
	L2= dlpfc.genes[["L2"]],
	L3.4= c(dlpfc.genes[["L3"]], dlpfc.genes[["L4"]]),
	L5= dlpfc.genes[["L5"]],
	L6= dlpfc.genes[["L6"]],
	WM= dlpfc.genes[["Oligo"]]
)

enrich.df_smooth$plot_genes = FALSE
for(i in names(sm_map)) {
  enrich.df_smooth$plot_genes = ifelse((enrich.df_smooth$gene_name %in% sm_map[[i]]) & 
                                     (enrich.df_smooth$smoothed==i), T, enrich.df_smooth$plot_genes)
}
## this is for the other version of the labeled volcano plot (see p3 in SZBDMulti-seq plot_layer-enrichment code
#enrich.df_smooth$gene_color = NA
#for(i in names(dlpfc.genes)) {
#  enrich.df_smooth$gene_color = ifelse((enrich.df_smooth$gene_name %in% dlpfc.genes[[i]]) & 
#                                     (enrich.df_smooth$plot_genes==T), i, enrich.df_smooth$gene_color)
#}

p4 <- ggplot(filter(enrich.df_smooth, gene_type_ptn==T), 
             aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point(data=filter(enrich.df_smooth, gene_type_ptn==T, plot_genes==F), color="grey70", size=.3)+
  geom_point(data=filter(enrich.df_smooth, gene_type_ptn==T, plot_genes==T), color="black", size=1)+
  facet_grid2(rows= vars(smoothed), scales="free_y",
              strip=strip)+
  ggrepel::geom_text_repel(data= filter(enrich.df_smooth, plot_genes==T),
                            aes(label=gene_name), 
                            size=3, fontface="bold.italic", 
                            min.segment.length=0, max.overlaps=Inf)+
  xlim(-5,5)+
  labs(title="PRECAST (smoothed)", subtitle="Protein coding genes only")+
  theme_bw()+theme(panel.grid.minor=element_blank(), strip.text.y=element_text(angle=0),
                   plot.margin = margin(.2,2.5,.2,2.5,"cm"))


pdf(file="plots/06_pseudobulk/PRECAST_smoothed/smoothed-k9-1663_layer-enrichment_plots.pdf",
    width=6, height=8)
p2
p4
p1
dev.off()

cat("\nPlots saved to: plots/06_pseudobulk/PRECAST_smoothed/smoothed-k9-1663_layer-enrichment_plots.pdf\n")


#top top layer markers
cat("\nTop layer markers (adj. p<1e-30, logFC>1\n")
filter(enrich.df_smooth, adj.P.Val<1e-30, sig_group=="large effect") %>%
  group_by(smoothed) %>% tally()


load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")

top.list = lapply(levels(enrich.df_smooth$smoothed), function(x) {
  tmp = filter(enrich.df_smooth, adj.P.Val<1e-30, sig_group=="large effect", smoothed==x)$gene_name
  phm = plotGroupedHeatmap(spe_pseudo, features=tmp, swap_rownames="gene_name",
                     group="smoothed_k9_1663", center=T, cluster_cols=F, angle_col=0, silent=T)
  phm[[4]]
  })
names(top.list) = levels(enrich.df_smooth$smoothed)

phm.list = gridExtra::marrangeGrob(top.list, ncol=1, nrow=1, top=quote(names(top.list)[g]))
ggsave("plots/06_pseudobulk/PRECAST_smoothed/smoothed-k9-1663_layer-markers-adjp-1e30-logfc-1_heatmap.pdf", phm.list, width=7, height=11)

cat("\n\nTop layer markers (adj p<1e-30, logFC>1) heatmap saved to: plots/06_pseudobulk/PRECAST_smoothed/smoothed-k9-1663_layer-markers-adjp-1e30-logfc-1_heatmap.pdf\n")


cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
