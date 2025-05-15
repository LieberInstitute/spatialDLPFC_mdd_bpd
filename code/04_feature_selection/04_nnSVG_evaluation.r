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
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")
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
             aes(x=factor(nsig_less10, levels=c("FALSE","TRUE"), labels=c("retain","Sig. SVG in\n<10 models")), 
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
             aes(x=factor(prop.sig_less33perc, levels=c("FALSE","TRUE"), labels=c("retain","Sig. SVG in\n<33% of models")), 
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
             aes(x=factor(best.rank_more500, levels=c("FALSE","TRUE"), labels=c("retain","Sig. SVG with\nbest rank >500")), 
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

geneList$qual_genes = include.genes


tmp.df = filter(results.df, gene_name %in% geneList$qual_genes) %>% mutate(is_qual= padj<.05 & rank<=500) %>%
  left_join(avg.expr)

p7 <- ggplot(left_join(tmp.df, filter(avg.expr, gene_name %in% geneList$qual_genes) %>% group_by(decile) %>% tally(name="n_genes_per_decile")),
       aes(x=rank, y=as.factor(decile), fill=n_genes_per_decile, lty=factor(is_qual, levels=c("TRUE","FALSE"))))+
  ggridges::geom_density_ridges(color="black")+
  scale_fill_gradientn(colors=c("grey90","lightblue","steelblue","navy"), values=c(0,.15,.4,.8,1))+
  labs(y="logcount expr decile\n(avg. over all 530k spots)", x="rank", 
       title=paste0("SVG feature list (n=",length(geneList$qual_genes),")"), fill="# top genes\nin decile",
       lty="qual. samples")+
  theme_minimal()+theme(axis.title.y=element_text(size=10))
#boxplot of spcov by decile of qual and nonqual genes
p8 <- ggplot(left_join(tmp.df, filter(avg.expr, gene_name %in% geneList$qual_genes) %>% group_by(decile) %>% tally(name="n_genes_per_decile")),
       aes(x=spcov, y=as.factor(decile), fill=n_genes_per_decile, lty=factor(is_qual, levels=c("TRUE","FALSE"))))+
  geom_boxplot(outliers=F)+
  scale_fill_gradientn(colors=c("grey90","lightblue","steelblue","navy"), values=c(0,.15,.4,.8,1))+
  labs(x="spatial coeff. of variance", y="logcount expr decile\n(avg. over all 530k spots)", 
       title=paste0("SVG feature list (n=",length(geneList$qual_genes),")"), fill="# top genes\nin decile",
       lty="qual. samples")+
  theme_minimal()+theme(axis.title.y=element_text(size=10))

lay.mat = rbind(c(1,1,2,2,3,3),c(4,4,4,5,5,5))
ggsave("plots/04_feature_selection/nnSVG-eval_find-candidate-SVGs.png", 
       gridExtra::grid.arrange(p3, p4, p5, #p6, 
		p7, p8, layout_matrix=lay.mat),
       bg="white", width=12, height=8, units="in")
cat("\nCandidate SVG filters plot(s) saved to: plots/04_feature_selection/nnSVG-eval_find-candidate-SVGs.png\n")

#layer marker heatmap
layer.markers = read.csv("processed-data/04_feature_selection/EXT_TableS9_sig_genes_FDR5perc_enrichment.csv") %>%
  filter(spatial_domain_resolution=="Sp09") %>%
  mutate(domain_simple=factor(test, 
                              levels=paste0("Sp09D0",c(1,2,3,5,8,4,7,6,9)), 
                              labels=c("L1 (1)","L1 (2)","L2","L3","L4","L5","L6","WM (1)","WM (2)")))

lm.svgs = filter(layer.markers, gene %in% geneList$qual_genes)
df1 = tidyr::pivot_wider(ungroup(lm.svgs)[,c("gene","stat","domain_simple")], 
                         names_from="domain_simple", values_from="stat",
                         values_fill = 0)
dim(df1)
colnames(df1)
m1 = as.matrix(df1[,-1])
rownames(m1) = df1$gene
hmp = pheatmap(m1[,c(1,2,3,5,8,4,7,6,9)], show_rownames = F, cluster_cols=F, angle_col = 0,
              main=paste(dim(m1)[[1]],"of 1663 SVGs are DLPFC layer marker genes\n(enrichment t stat fill color;\ngenes with FDR>.05 set to 0)"),
              )

ggsave(filename="plots/04_feature_selection/nnSVG-eval_1663-SVGs_layer-marker-heatmap.png",
       hmp[[4]], width=6, height=8.5, bg="white")
cat("\nLayer marker heatmap saved to: plots/04_feature_selection/nnSVG-eval_1663-SVGs_layer-marker-heatmap.png\n")

saveRDS(geneList, "processed-data/04_feature_selection/nnSVG-eval_geneList.rds")
length(unlist(geneList))
sapply(geneList, length)
cat("\nComplete list of all genes from 6111 that were filtered out has been saved to: processed-data/04_feature_selection/nnSVG-eval_geneList.rds\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
