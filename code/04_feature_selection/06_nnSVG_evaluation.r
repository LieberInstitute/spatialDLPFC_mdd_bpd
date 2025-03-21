setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(dplyr)
	library(ggplot2)
	library(pheatmap)
})
set.seed(123)

#load in supp files
avg.expr = read.csv("processed-data/04_feature_selection/nnSVG-filtered-genes_avg-logcounts.csv", row.names=1) %>%
  tibble::rownames_to_column(var="gene_id")

exclude.genes = readRDS("processed-data/04_feature_selection/batch-effect-genes_dummyslide-sample-seq-sex-condition_list.rds")
length(unique(unlist(exclude.genes))) #47
wm.keep = c("AQP1", "CERCAM", "NKX6-2", "MYRF")

exclude.names = filter(avg.expr, gene_id %in% unlist(exclude.genes))$gene_name
exclude.names2 = setdiff(exclude.names, wm.keep)

#load in nnSVG results
fileList = list.files("processed-data/04_feature_selection/per-sample_svgs")
length(fileList) #119!!!

resList = lapply(fileList, function(x) {
  name1 = substr(x, start=0, stop=13)
  df = read.csv(paste0("processed-data/04_feature_selection/per-sample_svgs/",x), row.names=1)
  rownames(df) <- NULL
  df$sample_id = name1
  return(df)
})
results.df = do.call(rbind, resList)

#switch to prop of spots
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n120_postQC_norm_")
results.df = left_join(results.df, as.data.frame(colData(spe)) %>% group_by(sample_id) %>% tally(name="n_spots_total")) %>%
  mutate(prop_spots_nonzero = n_spots_nonzero/n_spots_total)

#p value distribution and min number of spots per sample per gene -- NOT ADJUSTED P VAL
#http://varianceexplained.org/statistics/interpreting-pvalue-histogram/
pval.df = bind_rows(select(results.df, sample_id, gene_name, pval) %>% mutate(spots_nonzero=">100 spots"),
                    filter(results.df, prop_spots_nonzero>.05) %>% select(sample_id, gene_name, pval) %>% mutate(spots_nonzero="> 5%"),
                    filter(results.df, prop_spots_nonzero>.1) %>% select(sample_id, gene_name, pval) %>% mutate(spots_nonzero="> 10%")) %>%
  mutate(spots_nonzero=factor(spots_nonzero, levels=c(">100 spots","> 5%","> 10%"),
                                labels=c(">100 spots per sample",">5% of spots per sample",">10% of spots per sample")))

p1 <- ggplot(pval.df, aes(x=pval, group=sample_id))+
  stat_ecdf()+facet_grid(cols=vars(spots_nonzero))+
  scale_x_continuous(expand=expansion(add=c(0)))+
  theme_minimal()+theme(plot.margin=margin(0.5,1,.5,.5, unit="cm"),
                        panel.spacing = unit(1,"cm"),
                        strip.text=element_text(size=12), panel.ontop = T)

pval.df2 = left_join(pval.df, avg.expr) %>% group_by(spots_nonzero, decile, gene_name) %>% tally()
p2 <- ggplot(pval.df2, aes(x=as.factor(decile), y=n, color=as.factor(decile)))+
  ggbeeswarm::geom_quasirandom(size=.7)+facet_grid(cols=vars(spots_nonzero))+
  scale_color_viridis_d(option="turbo")+
  labs(x="logcount expr decile", y="# nnSVG models with gene")+
  theme_minimal()+theme(plot.margin=margin(0.5,1,.5,.5, unit="cm"),
                        panel.spacing = unit(1,"cm"),
                        legend.position="none", strip.text=element_text(size=12))

ggsave("plots/04_feature_selection/nnSVG-eval_pval-distribution_num-nonzero-spots.png", gridExtra::grid.arrange(p1, p2, ncol=1),
	bg="white", width=12, height=8, units="in")
cat("\nNon-zero spot filter plot saved to: plots/04_feature_selection/nnSVG-eval_pval-distribution_num-nonzero-spots.png\n")

#set # of nonzero spots and calculate FDR based on this new filtered set
results.df = filter(results.df, prop_spots_nonzero>.1) %>% group_by(sample_id) %>%
  mutate(padj=p.adjust(pval, method="fdr"))

geneList <- list()
cat("\nGenes removed because nonzero in <10% spots in every sample:\n")
geneList$less10perc_all = setdiff(avg.expr$gene_name, unique(results.df$gene_name))
length(geneList$less10perc_all)

cat("\nGenes removed because not sig in any nnSVG model:\n")
geneList$never_sig = group_by(results.df, gene_name) %>% summarise(n_run=n(), n_notsig=sum(padj>.05)) %>% 
  filter(n_run==n_notsig) %>% pull(gene_name)
length(geneList$never_sig)

#6106 starting 
cat("\nStarting gene set:\n")
include.genes = setdiff(avg.expr$gene_name, unique(unlist(geneList)))
length(include.genes)


cat("\nGenes sig in <10 nnSVG models:\n")
geneList$nsig_less10 = filter(results.df, gene_name %in% include.genes) %>%
  group_by(gene_name) %>% summarise(n_sig=sum(padj<.05)) %>% 
  filter(n_sig<10) %>% pull(gene_name)
length(geneList$nsig_less10)


p3 <- ggplot(filter(avg.expr, gene_name %in% include.genes) %>%
               mutate(nsig_less10=gene_name %in% geneList$nsig_less10), 
             aes(x=factor(nsig_less10, levels=c("FALSE","TRUE"), labels=c("retain","# sig models <10")), 
                 y=avg_expr))+
  ggbeeswarm::geom_quasirandom(width=.5, size=.5)+
  scale_y_continuous("logcount expr (avg. over all 530k spots)",
                     trans="log2", limits=c(0.005, 4), breaks=c(.03125, .25, 2.0),
                     labels=c("0.03125","0.25","2.0"))+
  labs(x="", title=paste("Remove",length(geneList$nsig_less10), "out of", length(include.genes),"remaining genes"))+
  theme_minimal()+theme(plot.title=element_text(size=11), plot.title.position="plot",
                        plot.margin=margin(.5,1,.5,.5, unit="cm"),
                        axis.text.x=element_text(size=12),
                        axis.title.y=element_text(size=10), axis.text.y=element_text(size=8))

#genes remaining
cat("\nGenes remaining:\n")
include.genes = setdiff(include.genes, unlist(geneList))
length(include.genes)

cat("Genes sig in <1/3 of models run:\n")
geneList$prop.sig_less33perc <- filter(results.df, gene_name %in% include.genes) %>% 
  group_by(gene_name) %>% summarise(n_run=n(), n_sig=sum(padj<.05)) %>%
  filter(n_sig/n_run<.33) %>% pull(gene_name)
length(geneList$prop.sig_less33perc)

p4 <- ggplot(filter(avg.expr, gene_name %in% include.genes) %>%
               mutate(prop.sig_less33perc=gene_name %in% geneList$prop.sig_less33perc), 
             aes(x=factor(prop.sig_less33perc, levels=c("FALSE","TRUE"), labels=c("retain","prop. sig models <0.33")), 
                 y=avg_expr))+
  ggbeeswarm::geom_quasirandom(width=.5, size=.5)+
  scale_y_continuous("logcount expr (avg. over all 530k spots)",
                     trans="log2", limits=c(0.005, 4), breaks=c(.03125, .25, 2.0),
                     labels=c("0.03125","0.25","2.0"))+
  labs(x="", title=paste("Remove",length(geneList$prop.sig_less33perc), "out of", length(include.genes),"remaining genes"))+
  theme_minimal()+theme(plot.title=element_text(size=11), plot.title.position="plot",
                        plot.margin=margin(.5,1,.5,.5, unit="cm"),
                        axis.text.x=element_text(size=12),
                        axis.title.y=element_text(size=10), axis.text.y=element_text(size=8))

#genes remaining
cat("\nGenes remaining:\n")
include.genes = setdiff(include.genes, unlist(geneList))
length(include.genes)


top.genes = filter(results.df, gene_name %in% include.genes, padj<.05, rank<=500) 
cat("\nGenes with best sig rank more than 500:\n")
geneList$best.rank_more500 = filter(results.df, gene_name %in% include.genes, !gene_name %in% unique(top.genes$gene_name)) %>%
  pull(gene_name) %>% unique()
length(geneList$best.rank_more500)


p5 <- ggplot(filter(avg.expr, gene_name %in% include.genes) %>%
               mutate(best.rank_more500=gene_name %in% geneList$best.rank_more500), 
             aes(x=factor(best.rank_more500, levels=c("FALSE","TRUE"), labels=c("retain","best sig. rank >500")), 
                 y=avg_expr))+
  ggbeeswarm::geom_quasirandom(width=.5, size=.5)+
  scale_y_continuous("logcount expr (avg. over all 530k spots)",
                     trans="log2", limits=c(0.005, 4), breaks=c(.03125, .25, 2.0),
                     labels=c("0.03125","0.25","2.0"))+
  labs(x="", title=paste("Remove",length(geneList$best.rank_more500), "out of", length(include.genes),"remaining genes"))+
  theme_minimal()+theme(plot.title=element_text(size=11), plot.title.position="plot",
                        plot.margin=margin(.5,1,.5,.5, unit="cm"),
                        axis.text.x=element_text(size=12),
                        axis.title.y=element_text(size=10), axis.text.y=element_text(size=8))

#1663 genes remaining
cat("\nGenes remaining:\n")
include.genes = setdiff(include.genes, unlist(geneList))
length(include.genes)


cat("\nNumber of total batch effect genes found:\n")
geneList$all_batch_effect = exclude.names2
length(geneList$all_batch_effect)
#cat("\nNumber of qualifying batch effect genes removed:\n")
#geneList$qual_batch_effect <- intersect(exclude.names2, include.genes)
#length(geneList$qual_batch_effect)
#sort(geneList$qual_batch_effect)

#p6 <- ggplot(filter(avg.expr, gene_name %in% include.genes) %>%
#         mutate(batch_effect=gene_name %in% geneList$qual_batch_effect),
#       aes(x=factor(batch_effect, levels=c("FALSE","TRUE"), labels=c("retain","batch effect")),
#           y=avg_expr))+
#  ggbeeswarm::geom_quasirandom(width=.5, size=.5)+
#  scale_y_continuous("logcount expr (avg. over all 530k spots)",
#                     trans="log2", limits=c(0.005, 4), breaks=c(.03125, .25, 2.0),
#                     labels=c("0.03125","0.25","2.0"))+
#  labs(x="", title=paste("Remove",length(geneList$qual_batch_effect), "out of", length(include.genes),"remaining genes"))+
#  theme_minimal()+theme(plot.title=element_text(size=11), plot.title.position="plot",
#                        plot.margin=margin(.5,1,.5,.5, unit="cm"),
#                        axis.text.x=element_text(size=12),
#                        axis.title.y=element_text(size=10), axis.text.y=element_text(size=8))

##1629 gene remaining
#cat("\nGenes remaining:\n")
#include.genes = setdiff(include.genes, unlist(geneList))
#length(include.genes)


geneList$qual_genes = include.genes


tmp.df = filter(results.df, gene_name %in% geneList$qual_genes) %>% mutate(is_qual= padj<.05 & rank<=500) %>%
  left_join(avg.expr)

p7 <- ggplot(left_join(tmp.df, filter(avg.expr, gene_name %in% geneList$qual_genes) %>% group_by(decile) %>% tally(name="n_genes_per_decile")),
       aes(x=rank, y=as.factor(decile), fill=n_genes_per_decile, lty=factor(is_qual, levels=c("TRUE","FALSE"))))+
  ggridges::geom_density_ridges(color="black")+
  scale_fill_gradientn(colors=c("grey90","lightblue","steelblue","navy"), values=c(0,.15,.4,.8,1))+
  labs(y="logcount expr decile\n(avg. over all 530k spots)", x="rank", 
       title=paste0("Candidate SVGs (n=",length(geneList$qual_genes),")"), fill="# top genes\nin decile",
       lty="qual. samples")+
  theme_minimal()+theme(axis.title.y=element_text(size=10))
#boxplot of spcov by decile of qual and nonqual genes
p8 <- ggplot(left_join(tmp.df, filter(avg.expr, gene_name %in% geneList$qual_genes) %>% group_by(decile) %>% tally(name="n_genes_per_decile")),
       aes(x=spcov, y=as.factor(decile), fill=n_genes_per_decile, lty=factor(is_qual, levels=c("TRUE","FALSE"))))+
  geom_boxplot(outliers=F)+
  scale_fill_gradientn(colors=c("grey90","lightblue","steelblue","navy"), values=c(0,.15,.4,.8,1))+
  labs(x="spatial coeff. of variance", y="logcount expr decile\n(avg. over all 530k spots)", 
       title=paste0("Candidate SVGs (n=",length(geneList$qual_genes),")"), fill="# top genes\nin decile",
       lty="qual. samples")+
  theme_minimal()+theme(axis.title.y=element_text(size=10))

lay.mat = rbind(c(1,1,2,2,3,3),c(4,4,4,5,5,5))
ggsave("plots/04_feature_selection/nnSVG-eval_find-candidate-SVGs.png", 
       gridExtra::grid.arrange(p3, p4, p5, #p6, 
		p7, p8, layout_matrix=lay.mat),
       bg="white", width=12, height=8, units="in")
cat("\nCandidate SVG filters plot(s) saved to: plots/04_feature_selection/nnSVG-eval_find-candidate-SVGs.png\n")

############ narrow down highly expressed (decile==10) genes
cat("\n\nNarrow down highly expressed (decile==10) genes\n")
geneList$top_decile = filter(avg.expr, gene_name %in% geneList$qual_genes, decile==10)$gene_name
cat("\nHighly expressed genes:\n")
length(geneList$top_decile)
length(geneList$top_decile)/length(geneList$qual_genes)

#load in correlation results
#correlations obtained from: code/04_feature_selection/06-supp_correlatePairs_decile10
corr_test = read.csv("processed-data/04_feature_selection/top-genes-padj05-rank500_decile10_correlation.csv")
colSums(is.na(corr_test))
#mirror for normal diagonal mtx
flip_test = corr_test
colnames(flip_test) = c("gene2","gene1","rho","p.value","FDR","gene2_name","gene2_avg_expr","gene1_name","gene1_avg_expr")
d1 = bind_rows(corr_test, flip_test)
length(unique(d1$gene1_name)) 

d2= tidyr::pivot_wider(d1[,c("gene1_name","gene2_name", "rho")], names_from="gene2_name", values_from="rho")
m1 = as.matrix(d2[,-1])
rownames(m1) = d2$gene1_name
#corr test was run on all highly expressed genes so make sure subset correctly
m2 = m1[geneList$top_decile, geneList$top_decile]

p9 <- pheatmap(m2, show_rownames=F, show_colnames=F, silent=T,
	treeheight_row=10, treeheight_col=10)

geneList$top.decile_high.spcov = filter(tmp.df, gene_name %in% geneList$top_decile) %>% group_by(gene_name) %>%
  summarise(med_spcov=median(spcov)) %>% filter(round(med_spcov,2)>=.4) %>%
  pull(gene_name)
geneList$top.decile_low.spcov = setdiff(geneList$top_decile, geneList$top.decile_high.spcov)

cat("\nKeeping highly expressed genes with spcov>=.40:\n")
length(geneList$top.decile_high.spcov)

p10 <- ggplot(filter(tmp.df, gene_name %in% geneList$top_decile) %>% 
               mutate(high.spcov=gene_name %in% geneList$top.decile_high.spcov), 
       aes(x=spcov, color=high.spcov, group=gene_name))+
  stat_ecdf()+scale_color_manual(values=c("black","red3"))+
  labs(x="spatial coefficient of variance (spcov)", title="Keep only high spcov (med>=0.40)")+
  theme_minimal()+theme(legend.position="inside", legend.position.inside = c(.85,.25),
                        legend.background = element_rect(fill="white", color="grey50"),
                        plot.margin=margin(.5,.5,0,.5, unit="cm"), plot.title.position="plot")

m3 = m2[geneList$top.decile_high.spcov,]
rownames(m3) = ifelse(rownames(m3) %in% geneList$all_batch_effect, paste(rownames(m3),"*"), rownames(m3))
p11 <- pheatmap(m3, show_rownames=T, show_colnames=F, 
               treeheight_row = 5, fontsize_row = 8, silent=T,
		treeheight_col=10)

ggsave("plots/04_feature_selection/nnSVG-eval_refine-candidate-SVGs_corr-ecdf.png", 
       gridExtra::grid.arrange(p9[[4]], p10, p11[[4]], ncol=3),
       bg="white", width=12, height=4, units="in")
cat("\nNarrowing highly expressed candidate SVG corr. plots saved to: plots/04_feature_selection/nnSVG-eval_refine-candidate-SVGs_corr-ecdf.png\n")

#dotplot 
#spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n120_postQC_norm_")
spe$slide2 = ifelse(spe$slide=="V13B23-283","V13B23-339",spe$slide)
spe$array2 = ifelse(spe$slide=="V13B23-283", "B_1", spe$array)
spe$sample_id2 = paste(spe$slide2, spe$array2)

spe_summ = scuttle::aggregateAcrossCells(spe[filter(avg.expr, gene_name %in% geneList$top.decile_high.spcov)$gene_id,],
                                         ids=spe$sample_id2, 
                                         statistics=c("mean","prop.detected"),
                                         use.assay.type="logcounts")

dotplot.df = left_join(tibble::rownames_to_column(as.data.frame(assay(spe_summ, "logcounts.mean")), var="gene_id") %>%
                         tidyr::pivot_longer(colnames(spe_summ), names_to="sample_id2", values_to="mean_expr"),
                       tibble::rownames_to_column(as.data.frame(t(scale(t(assay(spe_summ, "logcounts.mean"))))), var="gene_id") %>%
                         tidyr::pivot_longer(colnames(spe_summ), names_to="sample_id2", values_to="mean_expr_scaled")) %>%
  left_join(tibble::rownames_to_column(as.data.frame(assay(spe_summ, "logcounts.prop.detected")), var="gene_id") %>%
              tidyr::pivot_longer(colnames(spe_summ), names_to="sample_id2", values_to="prop_spots")) %>%
  left_join(avg.expr)
row.order = p11$tree_row$labels[p11$tree_row$order]
rev.row.order = rev(row.order) 
dotplot.df$y_order = factor(ifelse(dotplot.df$gene_name %in% geneList$all_batch_effect, 
	paste(dotplot.df$gene_name,"*"), dotplot.df$gene_name), levels=rev.row.order)

p12 = ggplot(dotplot.df, aes(x=sample_id2, y=y_order, color=mean_expr, size=prop_spots))+
  geom_count()+scale_color_viridis_c(option="F", direction=-1)+
  scale_y_discrete(position="right")+labs(size="Prop. spots")+
  scale_size(range=c(1,4), limits=c(0,1))+labs(color="Avg. expr.\n(logcounts)", y="")+
  theme_minimal()+theme(axis.text.x=element_blank(), legend.key.size = unit(14, "pt"),
                        plot.margin=margin(.5,.5,.5,.5, "cm"))

p13 = ggplot(dotplot.df, aes(x=sample_id2, y=y_order, color=mean_expr_scaled, size=prop_spots))+
  geom_count()+scale_color_viridis_c(option="F", direction=-1)+
  scale_y_discrete(position="right")+labs(size="Prop. spots")+
  scale_size(range=c(1,4), limits=c(0,1))+labs(color="Avg. expr.\n(scaled)", y="")+
  theme_minimal()+theme(axis.text.x=element_blank(), legend.key.size = unit(14, "pt"),
                        plot.margin=margin(.5,.5,.5,.5, "cm"))

ggsave("plots/04_feature_selection/nnSVG-eval_refine-candidate-SVGs_dotplot.png", 
       gridExtra::grid.arrange(p12, p13, ncol=1),
       bg="white", width=12, height=8, units="in")
cat("\nNarrowing highly expressed candidate SVG dot plots saved to: plots/04_feature_selection/nnSVG-eval_refine-candidate-SVGs_dotplot.png\n")

#remove any additional ribo genes
### without this extra step (before i fixed the FDR), there were 1053 bc of 3 RPS|RPL genes
### i found that the addition of those 3 genes really messed up precast (results saved in processed data)
#cat("\nGenes remaining after removing highly expressed genes with low spcov:\n")
#tmpList = setdiff(geneList$qual_genes, geneList$top.decile_low.spcov)
#length(tmpList)

#cat("\nRemaining RPS or RPL genes:\n")
#geneList$ribo <- grep("RPS|RPL", geneList$qual_genes, value=T)
#length(intersect(tmpList, geneList$ribo))

#final gene list
cat("\n\n#######################################################################################")
cat("\n################### FINAL SVG LIST\n\n")
geneList$final_svgs = setdiff(geneList$qual_genes, union(geneList$top.decile_low.spcov, geneList$all_batch_effect))
	#, geneList$ribo))
length(geneList$final_svgs)
cat("\n")
filter(avg.expr, gene_name %in% geneList$final_svgs) %>% group_by(decile) %>% tally()

p14 <- ggplot(filter(avg.expr, !gene_name %in% c(geneList$less10perc_all, geneList$never_sig)) %>%
         mutate(final.svgs=gene_name %in% geneList$final_svgs, cand.svgs=gene_name %in% geneList$qual_genes,
		final.svgs=ifelse(gene_name %in% geneList$all_batch_effect,"batch",as.character(final.svgs))), 
       aes(x=factor(final.svgs, levels=c("TRUE","FALSE","batch"), labels=c("retain","discard\n(other)","discard\n(batch)")), 
           y=avg_expr, color=cand.svgs))+
  ggbeeswarm::geom_quasirandom(width=.5, size=.5)+scale_color_manual("candidate\nSVGs",values=c("grey","black"))+
  scale_y_continuous("logcount expr (avg. over all 530k spots)",
                     trans="log2", limits=c(0.005, 4), breaks=c(.03125, .25, 2.0),
                     labels=c("0.03125","0.25","2.0"))+
  labs(x="", title=paste("Final SVG list:",length(geneList$final_svgs), "out of", nrow(avg.expr)-length(unlist(geneList[1:2])),"genes"))+
  theme_minimal()+theme(plot.title=element_text(size=11), plot.title.position="plot",
                        plot.margin=margin(.5,1,.5,.5, unit="cm"),
                        axis.text.x=element_text(size=12),
                        axis.title.y=element_text(size=10), axis.text.y=element_text(size=8))

#layer marker representation
layer.markers = read.csv("processed-data/04_feature_selection/EXT_TableS9_sig_genes_FDR5perc_enrichment.csv") %>%
  filter(stat>0, spatial_domain_resolution=="Sp09") %>%
  mutate(domain_simple=factor(test, 
                              levels=paste0("Sp09D0",c(1,2,3,5,8,4,7,6,9)), 
                              labels=c("L1","L1","L2","L3","L4","L5","L6","WM","WM")))
#collapse L1 and WM results to simple domain
layer.markers = group_by(layer.markers, domain_simple, gene) %>% slice_min(n=1, fdr) %>% ungroup()

final.svgs = filter(avg.expr, gene_name %in% geneList$final_svgs)
layerList = unique(layer.markers$domain_simple)
for(i in layerList) {
  final.svgs[[i]] = final.svgs$gene_name %in% filter(layer.markers, domain_simple==i)$gene
}

write.csv(final.svgs, paste0("processed-data/04_feature_selection/selected-SVGs_n",length(geneList$final_svgs),".csv"), row.names = F)
cat("\n\nFinal SVG dframe saved to:",paste0("processed-data/04_feature_selection/selected-SVGs_n",length(geneList$final_svgs),".csv"),"\n")

final.svgs2 = tidyr::pivot_longer(final.svgs, all_of(layerList), names_to="domain_simple", values_to="status") %>%
  filter(status==T) %>% mutate(domain_simple=as.character(domain_simple))

final.svgs2 = bind_rows(final.svgs2[,1:5], filter(final.svgs[,1:4], !gene_name %in% final.svgs2$gene_name) %>% mutate(domain_simple="none")) %>%
  mutate(domain_simple=factor(domain_simple, levels=c(levels(layerList),"none"))) %>%
  group_by(gene_name) %>% add_tally(name="n_sig_markers")

cat("\nRepresentation of layer markers in final SVG list:\n")
group_by(final.svgs2, domain_simple) %>% tally()

p15 <- ggplot(final.svgs2, aes(x=domain_simple, y=avg_expr, color=n_sig_markers))+
  ggbeeswarm::geom_quasirandom(width=.5, size=.5)+
  scale_color_gradientn("sig. marker for\nN domains",colors=c("#00204DFF","#6BAED6","#DD1C77"))+
  scale_y_continuous("logcount expr (avg. over all 530k spots)",
                     trans="log2", limits=c(0.005, 4), breaks=c(.03125, .25, 2.0),
                     labels=c("0.03125","0.25","2.0"))+
  labs(x="DLPFC layer marker", title=paste("Final SVG list:",length(geneList$final_svgs)))+
  theme_minimal()+theme(legend.key.size=unit(14, "pt"), plot.title=element_text(size=12), 
                        plot.margin=margin(.5,.5,.5,.5, unit="cm"),
                        axis.text.x=element_text(size=10),
                        axis.title.y=element_text(size=10), axis.text.y=element_text(size=8))


ggsave(paste0("plots/04_feature_selection/nnSVG-eval_final-SVG-list-n",length(geneList$final_svgs),"_dot-violin.png"),
       gridExtra::grid.arrange(p14, p15, layout_matrix=matrix(c(1,1,2,2,2), nrow=1)),
       bg="white", width=12, height=4, units="in")
cat("\nFinal SVG list dot/violin plots by layer markers saved to:",paste0("plots/04_feature_selection/nnSVG-eval_final-SVG-list-n",length(geneList$final_svgs),"_dot-violin.png"),"\n")


############ Look for highly expressed (decile==10) genes to supplement SVG list with
#cat("\n\n#######################################################################################")
#cat("\nLook for highly expressed (decile==10) genes to supplement SVG list with\n")
#
#qual_summary = filter(tmp.df, gene_name %in% geneList$top.decile_low.spcov) %>% group_by(gene_name) %>%
#  summarise(n_qual=sum(is_qual), n_sig=sum(padj<.05), n_notsig=sum(padj>.05),
#            n_sig_notqual=n_sig-n_qual)
#qual_summary$ribo = qual_summary$gene_name %in% grep("RPS|RPL", geneList$top.decile_low.spcov, value=T)
#
#p16 <- ggplot(qual_summary, aes(x=n_qual, y=n_notsig, color=ribo))+
#  geom_point()+geom_hline(aes(yintercept=9.5), color="red")+
#  scale_color_manual(values=c("grey","black"))+
#  labs(title=paste("Highly expressed genes with low spcov:",length(geneList$top.decile_low.spcov),"genes"),
#       x="# samples with gene ranked top 500 and padj<.05", y="# samples with gene padj>.05",
#       color="RPS | RPL")+
#  theme_minimal()+theme(plot.title=element_text(size=11), plot.title.position="plot", legend.title=element_text(size=10),
#                        plot.margin=margin(.5,.5,.5,.5, unit="cm"),
#                        axis.text=element_text(size=8),
#                        axis.title=element_text(size=10))
#
#cat("\nHighly expressed, low spcov genes with high # of not sig samples:\n")
#geneList$top.decile_low.spcov_high.notsig = filter(qual_summary, n_notsig>9)$gene_name
#length(geneList$top.decile_low.spcov_high.notsig)
#
#excluded.svgs = filter(avg.expr, gene_name %in% geneList$top.decile_low.spcov) %>% 
#  mutate(high_notsig=gene_name %in% geneList$top.decile_low.spcov_high.notsig)
#for(i in layerList) {
#  excluded.svgs[[i]] = excluded.svgs$gene_name %in% filter(layer.markers, domain_simple==i)$gene
#}
#
#excluded.svgs2 = tidyr::pivot_longer(excluded.svgs, all_of(layerList), names_to="domain_simple", values_to="status") %>%
#  filter(status==T) %>% mutate(domain_simple=as.character(domain_simple))
#
#excluded.svgs2 = bind_rows(excluded.svgs2[,1:6], filter(excluded.svgs[,1:5], !gene_name %in% excluded.svgs2$gene_name) %>% mutate(domain_simple="none")) %>%
#  mutate(domain_simple=factor(domain_simple, levels=c(levels(layerList),"none"))) %>%
#  group_by(gene_name) %>% add_tally(name="n_sig_markers")
#
#p17 <- ggplot(excluded.svgs2, aes(x=domain_simple, y=avg_expr, color=high_notsig))+
#  ggbeeswarm::geom_quasirandom(width=.5, size=.5)+scale_color_manual(values=c("black","red"))+
#  scale_y_continuous("logcount expr (avg. over all 530k spots)",
#                     trans="log2", limits=c(0.005, 4), breaks=c(.03125, .25, 2.0),
#                     labels=c("0.03125","0.25","2.0"))+
#  labs(x="DLPFC layer marker", color=">=10 samples\nwhere padj>.05",
#       title=paste("Highly expressed genes with low spcov:",length(geneList$top.decile_low.spcov),"genes"))+
#  theme_minimal()+theme(plot.title=element_text(size=12), 
#                        plot.margin=margin(.5,.5,.5,.5, unit="cm"),
#                        axis.text.x=element_text(size=10),
#                        axis.title.y=element_text(size=10), axis.text.y=element_text(size=8))
#
#ggsave("plots/04_feature_selection/nnSVG-eval_supplemental-gene-options.png",
#       gridExtra::grid.arrange(p16, p17, layout_matrix=matrix(c(1,2,2), nrow=1)),
#       bg="white", width=12, height=4, units="in")
#cat("\nHighly expressed supplemental gene plots saved to: plots/04_feature_selection/nnSVG-eval_supplemental-gene-options.png\n\n")


saveRDS(geneList, "processed-data/04_feature_selection/nnSVG-eval_geneList.rds")
length(unlist(geneList))
sapply(geneList, length)
cat("\nComplete list of all genes from 6111 that were filtered out has been saved to: processed-data/04_feature_selection/nnSVG-eval_geneList.rds\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
