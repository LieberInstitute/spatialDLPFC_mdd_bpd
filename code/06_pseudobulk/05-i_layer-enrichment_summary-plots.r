library(SpatialExperiment)
library(dplyr)
library(ggplot2)

load("processed-data/06_pseudobulk/spe_n119_pseudo_sample-comb-clus_norm-filt.Rdata")

avg.expr = read.csv("processed-data/06_pseudobulk/filtered-genes_avg-logcounts.csv", row.names=1)
layer_res = read.csv("processed-data/06_pseudobulk/results_layer-enrichment_covars-age-detected-ncells-sex-slide.csv", row.names=1)

#format layer enrichment results
tstats.df = layer_res[,c(grep("t_stat",colnames(layer_res)),33:34)] %>% 
  tidyr::pivot_longer(cols=grep("t_stat",colnames(layer_res), value=T), names_to="combined_cluster", values_to="t_stat", names_prefix="t_stat_")
fdr.df = layer_res[,c(grep("fdr",colnames(layer_res)),33:34)] %>% 
  tidyr::pivot_longer(cols=grep("fdr",colnames(layer_res), value=T), names_to="combined_cluster", values_to="fdr", names_prefix="fdr_")
logFC.df = layer_res[,c(grep("logFC",colnames(layer_res)),33:34)] %>% 
  tidyr::pivot_longer(cols=grep("logFC",colnames(layer_res), value=T), names_to="combined_cluster", values_to="logFC", names_prefix="logFC_")

enrich.df = left_join(tstats.df, fdr.df, by=c("combined_cluster","ensembl","gene")) %>%
  left_join(logFC.df, by=c("combined_cluster","ensembl","gene")) %>%
  mutate(combined_cluster=factor(combined_cluster, levels=c("Vasc","L1","L2","L3","GABA","L5","L6","WM")))


enrich.df$fdr_bin = cut(enrich.df$fdr, breaks=c(0,.0001,.05,1))
table(enrich.df$fdr_bin, useNA="ifany")
enrich.df$fdr_bin2 = as.character(factor(enrich.df$fdr_bin, levels=levels(enrich.df$fdr_bin), labels=c("highly sig.","sig.","NS")))
enrich.df$sig_group = ifelse(enrich.df$logFC>1, "large effect", enrich.df$fdr_bin2)
table(enrich.df[,c("sig_group","fdr_bin2")], useNA="ifany")
enrich.df$sig_group2 = ifelse(abs(enrich.df$logFC)>1, "large effect", enrich.df$fdr_bin2)
table(enrich.df[,c("sig_group2","fdr_bin2")], useNA="ifany")

#volcano plot
ggplot(enrich.df, aes(x=logFC, y=-log10(fdr), color=sig_group2))+
  geom_point(size=.5)+facet_wrap(vars(combined_cluster), ncol=2)+
  xlim(-6.2,6.2)+
  scale_color_manual(values=c("NS"="grey","sig."="grey50","highly sig."="black", "large effect"="red3"))+
  theme_bw()

#bar plot of large effect genes
tmp = filter(enrich.df, sig_group2=="large effect") %>% mutate(direction= sign(logFC)) %>% 
  group_by(combined_cluster, direction) %>% tally()
ggplot(tmp, aes(x=combined_cluster, y=n))+
  geom_bar(data=filter(tmp, direction==1), stat="identity", fill="tomato3")+
  geom_bar(data=filter(tmp, direction==-1) %>% mutate(n=-n), stat="identity", fill="skyblue")+
  labs(y="# sig. markers\nFDR<.0001, abs(logFC)>1", x="")+
  theme_bw()+theme(aspect.ratio=.8, text=element_text(size=14))

domains = levels(enrich.df$combined_cluster)
names(domains) = domains
upList = lapply(domains, function(x) {
  filter(enrich.df, combined_cluster==x, fdr<.0001, logFC>0)$gene
})
sapply(upList, length)

#overlap list identifies layer-specific markers (fdr<.0001, logFC>1) that are not present in other clusters (fdr<.0001, logFC>0)
overlapList =  lapply(domains, function(x) {
  query1 = filter(enrich.df, combined_cluster==x, gene %in% upList[[x]], logFC>1)$gene
  others = upList[setdiff(domains, x)]
  out1 = lapply(others, function(y) intersect(y, query1))
  out1$unique = setdiff(query1, unlist(others))
  return(out1)
})
str(overlapList)

markerList = lapply(overlapList, function(x) x$unique)
sapply(markerList, length)

#used this filter and scater expression plots to explore top markers to highlight in figure panels
filter(enrich.df, gene %in% markerList[["WM"]], combined_cluster=="WM") %>% 
  arrange(desc(t_stat)) %>% select(gene, combined_cluster, t_stat, fdr, logFC) %>%
  left_join(avg.expr, by=c("gene"="gene_name")) #%>%
  filter(gene=="TLX3")

scater::plotExpression(spe_pseudo, features=c("PCP4", "HTR2C","TRABD2A","TOX"), 
                       x="combined_cluster", swap_rownames = "gene_name", assay.type = "logcounts")

#marker genes to highlight in heatmap
marker.genes = c("MYL9","PDGFRB","CLDN5",
                 "FABP7","RELN","NCAN",#or AQP4 for astrocytes
                 "ITPKA","DGKA",#L2 only
		 "HPCAL1","RASGRF2",#both L2/3 that are fdr<.0001 in L2 only but also sig in L3
                 "KCNH5","TNNT2",#L3 only
                 "RORB",
                 "NXPH1","SLC32A1","ARX",#GABA
                 "PCP4", "HTR2C","TRABD2A","TOX",
		 "CLSTN2",
                 "MOXD1","ISLR2","OPRK1",
                 "SHTN1","CNP", "SLC44A1"
                 )

hmp = scater::plotGroupedHeatmap(spe_pseudo, features=marker.genes, group="combined_cluster", 
                           swap_rownames="gene_name", scale = T, center=T, cluster_cols=F, cluster_rows=F,
                           angle_col=0)

ggsave("plots/06_pseudobulk/test_layer-enrichment_heatmap.png",
       hmp[[4]], bg="white", width=4, height=6)

#subset to highlight in boxplots
marker.genes2 = c("MYL9",#Vasc
  "FABP7",
  "ITPKA",
  "KCNH5",
  "SLC32A1",
  "TOX",#"HTR2C",
  "ISLR2",#"OPRK1",
  "SHTN1" #CNP IS THE THIRD MOST ABUNDANT PROTEIN IN CNS MYELIN
  )

cdata = as.data.frame(colData(spe_pseudo))
for(i in marker.genes2) {
  cdata[[i]] = logcounts(spe_pseudo)[rowData(spe_pseudo)$gene_name==i,]
}

p1 <- ggplot(tidyr::pivot_longer(cdata[,c("sample_id","combined_cluster",marker.genes2)], all_of(marker.genes2),
                           names_to="gene", values_to="expr") %>%
         mutate(gene=factor(gene, levels=marker.genes2)), 
       aes(x=combined_cluster, y=expr, fill=combined_cluster))+
  geom_boxplot(outliers=F, linewidth=.3, width=.5)+#geom_violin(trim=F)+
  facet_wrap(vars(gene), ncol=1, scales="free_y")+
  scale_fill_manual(values=precast.colorList[["n1663_k9"]]$colors[1:8])+
  labs(y="log2 CPM")+
  theme_bw()+theme(legend.position="none", axis.title.x=element_blank(),
    strip.background = element_rect(fill="transparent", color="transparent"), 
    text=element_text(size=14), strip.text=element_text(size=10, face="italic"))

ggsave("plots/06_pseudobulk/test_layer-enrichment_boxplot.png",
       p1, bg="white", width=3.5, height=7)


p2 <- ggplot(tidyr::pivot_longer(cdata[,c("sample_id","combined_cluster",marker.genes2)], all_of(marker.genes2),
                                 names_to="gene", values_to="expr") %>%
               mutate(gene=factor(gene, levels=marker.genes2)), 
             aes(x=combined_cluster, y=expr, color=combined_cluster))+
  ggbeeswarm::geom_quasirandom(size=.5)+
  geom_boxplot(color="black", linewidth=.5, fill="transparent", outliers=F, width=.6)+
  facet_wrap(vars(gene), ncol=1, scales="free_y")+
  scale_color_manual(values=precast.colorList[["n1663_k9"]]$colors[1:8])+
  labs(y="log2 CPM")+
  theme_bw()+theme(legend.position="none", axis.title.x=element_blank(),
                   strip.background = element_rect(fill="transparent", color="transparent"), 
                   text=element_text(size=10), strip.text=element_text(size=10, face="italic"))

ggsave("plots/06_pseudobulk/test_layer-enrichment_violin.png",
       p2, bg="white", width=3.5, height=7)

#spot plots to highlight
library(HDF5Array)
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n120_postQC_norm_")
clusters = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
identical(rownames(colData(spe)), rownames(clusters))
spe$combined_cluster= factor(clusters$combined_cluster, levels=c("Vasc","L1","L2","L3","GABA","L5","L6","WM"))
spe$is_gaba = spe$combined_cluster=="GABA"
#(n16.samples = paste0("V13B23-30",c(1,2,8,9)))
#spe_sub = spe[,spe$slide %in% n16.samples]
spe_sub = spe[,spe$sample_id=="V13B23-308_D1" & !is.na(spe$combined_cluster)]

spot.genes = c("PCDH8","RORB","PCP4","CLSTN2",#general gradients
               "DGKA","TNNT2","TRABD2A","OPRK1")#specific 
for(i in spot.genes) {
  spe_sub[[i]] = logcounts(spe_sub)[rowData(spe_sub)$gene_name==i,]
}

plist <- lapply(spot.genes, function(x) {
  p = make_escheR(spe_sub) |> #add_ground(var="combined_cluster") |> 
    add_fill(var=x, point_size = 1)
  p+#scale_color_manual("", values=precast.colorList[["n1663_k9"]][["colors"]][1:8])+
    scale_fill_gradient("log2\nCPM", low="white",high="black")+
    labs(title=x)+theme(text=element_text(size=10), plot.title=element_text(face="italic"),
                        legend.position="none")
})

ggsave("plots/06_pseudobulk/test_layer-enrichment_spot-plots.png", 
       gridExtra::grid.arrange(plist[[1]], plist[[2]], plist[[3]], plist[[4]], 
                               plist[[5]], plist[[6]], plist[[7]], plist[[8]],
                               layout_matrix=cbind(1:4,5:8)),
       bg="white", width=6, height=8)


library(escheR)
source("code/05_clustering/PRECAST/PRECAST_colorLists.r")

