setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(DelayedArray)
	library(spatialLIBD)
	library(dplyr)
	library(ggplot2)
	library(pheatmap)
})

set.seed(123)

load("processed-data/06_pseudobulk/spe_n119_pseudo_sample-comb-clus_norm-filt.Rdata")
dim(spe_pseudo) # 12397   942

#layer_mod <- registration_model(spe_pseudo,
#	covars = c("detected","sex","slide"),
#	var_registration = "combined_cluster"
#)
#
#layer_block_cor <- registration_block_cor(spe_pseudo, registration_model = layer_mod,
#    var_sample_id = "sample_id"
#)
#
#layer_res <- registration_stats_enrichment(spe_pseudo, block_cor = layer_block_cor,
#  covars = c("detected","sex","slide"),
#  var_registration = "combined_cluster",
#  gene_ensembl = "gene_id",
#  gene_name = "gene_name"
#)

#write.csv(layer_res, "processed-data/06_pseudobulk/results_layer-enrichment_covars-detected-sex-slide.csv", row.names=T)
#cat("\nLayer enrichment test results saved to: processed-data/06_pseudobulk/results_layer-enrichment_covars-detected-sex-slide.csv\n")

################################################### actual spatial registration
layer_res = read.csv("processed-data/06_pseudobulk/results_layer-enrichment_covars-detected-ncells-sex-slide.csv", row.names=1)

#follow spatialLIBD spatial registration tutorial 
layer_modeling_results <- fetch_data(type = "modeling_results")

query_t_stats <- layer_res[, grep("^t_stat", colnames(layer_res))]
colnames(query_t_stats) <- gsub("^t_stat_", "", colnames(query_t_stats))

cor_layer <- layer_stat_cor(
  stats = query_t_stats,
  modeling_results = layer_modeling_results,
  model_type = "enrichment",
  top_n = 100
)

ggsave("plots/06_pseudobulk/layer-enrich_spatial-registration_covars-detected-ncells-sex-slide.png",
	layer_stat_cor_plot(cor_layer[c("Vasc","L1","L2","L3","GABA","L5","L6","WM"),], max = max(cor_layer)),
	bg="white", height=7, width=7, units="in"
)
cat("\nSpatial registration plot saved to: plots/06_pseudobulk/layer-enrich_spatial-registration_covars-detected-ncells-sex-slide.png\n")
################################################### 

#create cluster marker csv and explore top markers
source("code/05_clustering/PRECAST/PRECAST_colorLists.r")
avg.expr = read.csv("processed-data/06_pseudobulk/filtered-genes_avg-logcounts.csv", row.names=1)

tstats.df = layer_res[,c(grep("t_stat",colnames(layer_res)),33:34)] %>% 
  tidyr::pivot_longer(cols=grep("t_stat",colnames(layer_res), value=T), names_to="combined_cluster", values_to="t_stat", names_prefix="t_stat_")
fdr.df = layer_res[,c(grep("fdr",colnames(layer_res)),33:34)] %>% 
  tidyr::pivot_longer(cols=grep("fdr",colnames(layer_res), value=T), names_to="combined_cluster", values_to="fdr", names_prefix="fdr_")
logFC.df = layer_res[,c(grep("logFC",colnames(layer_res)),33:34)] %>% 
  tidyr::pivot_longer(cols=grep("logFC",colnames(layer_res), value=T), names_to="combined_cluster", values_to="logFC", names_prefix="logFC_")

enrich.df = left_join(tstats.df, fdr.df, by=c("combined_cluster","ensembl","gene")) %>%
  left_join(logFC.df, by=c("combined_cluster","ensembl","gene")) %>%
  mutate(combined_cluster=factor(combined_cluster, levels=c("Vasc","L1","L2","L3","GABA","L5","L6","WM")))

cat("\nNumber of layer markers with fdr<.0001, t_stat>0\n")
filter(enrich.df, fdr<.0001, t_stat>0) %>% group_by(combined_cluster) %>% tally()

cat("\nNumber of layer markers with fdr<.0001, t_stat>0, logFC>1\n")
filter(enrich.df, fdr<.0001, t_stat>0, logFC>1) %>% group_by(combined_cluster) %>% tally()

cat("\nHighly sig (fdr<.0001) and high logFC (>1) layer markers are spread across expression decile\n")
filter(enrich.df, fdr<.0001, t_stat>0) %>% mutate(logfc_more1=logFC>1) %>% group_by(gene) %>%
  summarise(n_clus_logfc_more1 = sum(logfc_more1)) %>% left_join(avg.expr, by=c("gene"="gene_name")) %>%
  select(n_clus_logfc_more1, decile) %>% table()

source.df = filter(enrich.df, fdr<.0001, t_stat>0, logFC>1) %>% group_by(gene) %>% add_tally(name="n_clus_sig") 
clus1.df = filter(source.df, n_clus_sig==1) %>% transmute(sig_clus=combined_cluster)
clus2.df = filter(source.df, n_clus_sig==2) %>% summarise(sig_clus=paste(combined_cluster, collapse="/"))
cluster.markers = bind_rows(clus1.df, filter(clus2.df, sig_clus %in% c("L1/L2", "L1/WM", "L2/GABA", "L2/L3", "L3/GABA", "L6/WM", "Vasc/L1", "Vasc/WM"))) %>%
	#c("L1/WM", "L2/L3", "L2/GABA","L6/WM", "Vasc/L1", "Vasc/WM"))) %>%
  left_join(avg.expr, by=c("gene"="gene_name")) %>%
  mutate(sig_clus=factor(sig_clus, levels=c("Vasc","Vasc/L1", "L1", "L1/L2", "L2", "L2/L3", "L3", "L5", "L6", "L6/WM", "WM", "L1/WM", "Vasc/WM", "L2/GABA", "L3/GABA", "GABA")))
	#c("Vasc","Vasc/L1","L1","L2","L2/GABA","L2/L3","L3","GABA","L5","L6","L6/WM","WM","Vasc/WM","L1/WM")))
table(cluster.markers$sig_clus)

write.csv(cluster.markers, "processed-data/06_pseudobulk/layer-enrich_covars-detected-ncells-sex-slide_fdr-0001-logfc-1_cluster-markers.csv", row.names=F)
cat("\nCluster markers dataframe saved to: processed-data/06_pseudobulk/layer-enrich_covars-detected-ncells-sex-slide_fdr-0001-logfc-1_cluster-markers.csv\n")

#plot top markers
plot.genes = c(filter(cluster.markers, !sig_clus %in%  c("Vasc","L1","GABA","WM"))$gene,
               filter(source.df, n_clus_sig==1, combined_cluster %in%  c("Vasc","L1","GABA","WM")) %>% group_by(combined_cluster) %>% slice_max(n=40, t_stat) %>% pull(gene)
          )
cat("\nNumber of top markers to plot:\n")
length(plot.genes)
plot.ids = rownames(spe_pseudo)[rowData(spe_pseudo)$gene_name %in% plot.genes]
row_annot = left_join(cbind.data.frame("gene_id"=plot.ids,"gene_name"=rowData(spe_pseudo)[plot.ids,"gene_name"]), cluster.markers, by=c("gene_name"="gene")) %>%
  arrange(sig_clus)
row_annot_df = data.frame("sig_clus"=as.character(row_annot$sig_clus))
rownames(row_annot_df) = row_annot$gene_id
annot_colors=list("sig_clus"=c(precast.colorList[["n1663_k9"]][["colors"]][1:8],
                               "#FDBF6F", "lightgreen", "#008080", "#8B0000", "#EF5327", "navy","#6A3D9A","#4B006E"))
names(annot_colors$sig_clus) = c("Vasc","L1", "L2", "L3",  "GABA", "L5", "L6", "WM",
                 "Vasc/L1","L1/L2", "L2/L3","L6/WM", "Vasc/WM", "L1/WM", "L2/GABA", "L3/GABA")
#annot_colors=list("sig_clus"=c(precast.colorList[["n1663_k9"]][["colors"]][1:8],
#                               "#FDBF6F","#6A3D9A","#008080", "#8B0000", "#EF5327", "navy"))
#names(annot_colors$sig_clus) = c("Vasc","L1","L2","L3","GABA","L5","L6","WM",
#                                 "Vasc/L1","L2/GABA","L2/L3","L6/WM","Vasc/WM","L1/WM")
annot_colors$sig_clus = annot_colors$sig_clus[levels(row_annot$sig_clus)]
groupList = list(ntc.f = spe_pseudo$condition=="NTC" & spe_pseudo$sex=="F",
     ntc.m = spe_pseudo$condition=="NTC" & spe_pseudo$sex=="M",
     mdd.f = spe_pseudo$condition=="MDD" & spe_pseudo$sex=="F",
     mdd.m = spe_pseudo$condition=="MDD" & spe_pseudo$sex=="M",
     bpd.f = spe_pseudo$condition=="BPD" & spe_pseudo$sex=="F",
     bpd.m = spe_pseudo$condition=="BPD" & spe_pseudo$sex=="M")

m1 = do.call(cbind, lapply(names(groupList), function(x) {
  clusList = c("Vasc","L1","L2","L3","GABA","L5","L6","WM")
  names(clusList) = clusList
  outList = lapply(clusList, function(y) {
    x1 = groupList[[x]] & spe_pseudo$combined_cluster==y
    rowMeans(logcounts(spe_pseudo)[plot.ids,x1])
  })
  names(outList) = paste(x, names(outList), sep="_")
  do.call(cbind, outList)
}))

col_annot = data.frame("condition"=substr(colnames(m1), start=0, stop=3))
rownames(col_annot) = colnames(m1)
annot_colors$condition = c(ntc="black",mdd="#9e771b", bpd="#1b9e77")
hmp = pheatmap(m1[row_annot$gene_id,], cluster_rows = F, 
                   annotation_row = row_annot_df, annotation_col= col_annot, annotation_colors = annot_colors,
                   scale="row", show_rownames = F)

ggsave("plots/06_pseudobulk/layer-enrich_top-cluster-markers_covars-detected-ncells-sex-slide.png", hmp[[4]], bg="white", height=7, width=7, units="in")
cat("\nTop layer markers heatmap saved to: plots/06_pseudobulk/layer-enrich_top-cluster-markers_covars-detected-ncells-sex-slide.png\n")


cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
