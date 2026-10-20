setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SingleCellExperiment)
	library(edgeR)
	library(dplyr)
	library(ggplot2)
})

set.seed(123)
source("code/07_dx_DE/custom_functions.r")
cpList = readRDS("plots/colorPalettes.rds")

results_set = "Maitra (MDD)"
load("processed-data/08_compare_DEGs/Maitra/sce_Maitra_pseudo_sample-low-res_norm-filt.Rdata")
spe_pseudo <- sce_pseudo # for use of custom functions

#comp_names = c("MicrodotVasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo")
#names(comp_names) = c("Micro.Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo")

comp_names = c("MdotV","Ast","L2","L3","L4","Inb","L5","L6","Olg")
names(comp_names) = c("M.V","Ast","L2","L3","L4","Inb","L5","L6","Olg")


comparisons = c("F_NTC.MDD","M_NTC.MDD")

### layer-adjusted 
results <- readRDS("processed-data/08_compare_DEGs/Maitra/lmFit-voom_layer-adjusted_low-res_condition-sex_covars-pc3-ncells.rda")

eb <- eBayes(groupContrasts(results, comparisons), trend=T, robust=T)
#plotSA(eb)

## F statistics
eb.f = getTopTable(eb, .coef="all") %>%
	mutate("F_design"="~ 0 + dx*sex")
f_sig = filter(eb.f, adj.P.Val<.05)$gene_id
## moderated t statistics
sex.res = sexTopTable(eb, phist=T)
### implement global padj method (without this, different t cutoffs were sig because of different group sizes)
### details: https://github.com/cran/limma/blob/master/R/decidetests.R
sex.res$results$adj.P.Val = p.adjust(sex.res$results$P.Value, method="BH")

cat("\nLayer-adjusted analysis with adj p<.05 and F adj p<.05:\n\n")
filter(sex.res$results, adj.P.Val<.05, gene_id %in% f_sig) %>% 
	group_by(sex, group, .drop=F) %>% tally()

p <- ggplot(sex.res$results, aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point(size=.3)+
  geom_point(data=filter(sex.res$results, adj.P.Val<.05, gene_id %in% f_sig), size=.3, color="red2")+
  facet_grid(rows=vars(group), cols=vars(sex))+
  xlim(-ceiling(max(sex.res$results$logFC)), ceiling(max(sex.res$results$logFC)))+
  ggtitle(paste0("Layer-adjusted (", results_set, ")"))+
  theme_bw()


ggsave("plots/08_compare_DEGs/Maitra/layer-adjusted-pc3-ncells_low-res_p-val-histogram-volcano.png",
	gridExtra::grid.arrange(sex.res$phist, p, ncol=1),
	bg="white", width=6, height=6)
cat("\n\nSaved un-adjusted p value histogram and volcano plots to: plots/08_compare_DEGs/Maitra/layer-adjusted-pc3-ncells_low-res_p-val-histogram-volcano.png\n")

write.csv(eb.f, "processed-data/08_compare_DEGs/Maitra/layer-adjusted-pc3-ncells_low-res_F-test.csv", row.names=F)
cat("\n\n\nSaved F test results dframe to: processed-data/08_compare_DEGs/Maitra/layer-adjusted-pc3-ncells_low-res_F-test.csv\n")
write.csv(sex.res$results, "processed-data/08_compare_DEGs/Maitra/layer-adjusted-pc3-ncells_low-res_moderated-t-test.csv", row.names=F)
cat("\nSaved moderated t test results dframe to: processed-data/08_compare_DEGs/Maitra/layer-adjusted-pc3-ncells_low-res_moderated-t-test.csv\n\n")


### layer restricted
results <- readRDS("processed-data/08_compare_DEGs/Maitra/lmFit-voom_layer-restricted_low-res_condition-sex_covars-pc3-ncells.rda")

#need to have "dot" instead of "." in comparisons for cluster name 
comparisons = unlist(lapply(comp_names, 
                            function(x) paste(x, comparisons, sep="_")))

eb <- eBayes(groupContrasts(results, comparisons), trend=T, robust=T)
#plotSA(eb)

## F statistics
eb_f = getTopTable(eb, .coef="all") %>%
	mutate(cluster="all", F_design="~ 0 + dx*sex*cluster")

## moderated t statistics
sex.resList <- lapply(names(comp_names), function(x) sexTopTable(eb, phist=F, cluster=x))
sex.res_df = do.call(rbind, sex.resList) %>%
  mutate(cluster=factor(cluster, levels=names(comp_names)))
### implement global padj method (without this, different t cutoffs were sig because of different group sizes)
### details: https://github.com/cran/limma/blob/master/R/decidetests.R
sex.res_df$adj.P.Val = p.adjust(sex.res_df$P.Value, method="BH")

#alt pvalue histogram
tmp = mutate(sex.res_df, sex.group=factor(paste(group, sex), levels=c("NTC.MDD F","NTC.MDD M")))
p0 <- ggplot(tmp, aes(x=P.Value))+
  geom_histogram(bins=30)+scale_x_continuous(breaks=c(0,.5,1))+
  facet_grid(rows=vars(cluster), cols=vars(sex.group), scales="free_y")+
  labs(title=paste0("Layer-restricted (", results_set, ")"))+
  theme_minimal()+theme(panel.grid.minor=element_blank())


#number of DEGs
cat("\nLayer-restricted analysis with adj p<.05 and F adj p<.05:\n\n")
filter(sex.res_df, adj.P.Val<.05, gene_id %in% eb_f$gene_id[eb_f$adj.P.Val<.05]) %>% 
  group_by(sex, group, cluster) %>% tally() %>%
  tidyr::pivot_wider(names_from="cluster", values_from="n", values_fill=0)

#volcanos
tmp = do.call(rbind, lapply(names(comp_names), function(x)
  filter(sex.res_df, cluster==x, adj.P.Val<.05, gene_id %in% eb_f$gene_id[eb_f$adj.P.Val<.05])
))

p1 <- ggplot(mutate(sex.res_df, facets=factor(paste(sex, group, sep="_"), levels=c("F_NTC.MDD","M_NTC.MDD"))), 
    aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point(size=.3)+geom_point(data=filter(tmp, sex=="F"), size=.3, color="red2")+
  facet_grid(rows=vars(cluster), cols=vars(facets))+
  xlim(-ceiling(max(sex.res_df$logFC)), ceiling(max(sex.res_df$logFC)))+ggtitle("Maitra 2023: MDD")+
  theme_bw()


pdf(file="plots/08_compare_DEGs/Maitra/layer-restricted-pc3-ncells_low-res_p-val-histogram-volcano.pdf",
    width=5, height=8)
p0
ggrastr::rasterize(p1,layer='point',dpi=300)
dev.off()
cat("\n\nSaved volcano plots to: plots/08_compare_DEGs/Maitra/layer-restricted-pc3-ncells_low-res_p-val-histogram-volcano.pdf\n\n\n")


write.csv(eb_f, "processed-data/08_compare_DEGs/Maitra/layer-restricted-pc3-ncells_low-res_F-test_all.csv", row.names=F)
cat("\n\n\nSaved F test results dframe to: processed-data/08_compare_DEGs/Maitra/layer-restricted-pc3-ncells_low-res_F-test_all.csv\n")
write.csv(sex.res_df, "processed-data/08_compare_DEGs/Maitra/layer-restricted-pc3-ncells_low-res_moderated-t-test.csv", row.names=F)
cat("\nSaved moderated t test results dframe to: processed-data/08_compare_DEGs/Maitra/layer-restricted-pc3-ncells_low-res_moderated-t-test.csv\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()

