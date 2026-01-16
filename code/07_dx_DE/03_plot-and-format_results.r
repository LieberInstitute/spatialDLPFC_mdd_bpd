setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(edgeR)
	library(dplyr)
	library(ggplot2)
})

set.seed(123)
source("code/07_dx_DE/custom_functions.r")
cpList = readRDS("plots/colorPalettes.rds")

#load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
#results_set = "smoothed-k9-1663"
#comp_names = c("L1","L2","L3dot4","L5","L6","WM")
#names(comp_names) = c("L1","L2","L3.4","L5","L6","WM")

load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
results_set = "seurat-pc30"
comp_names = c("MicrodotVasc","Astro","L2dot3","L4","Inhb","L5","L6","Oligo")
names(comp_names) = c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo")

## layer adjusted
results <- readRDS(paste0("processed-data/07_dx_DE/lmFit-voom_layer-adjusted_", results_set, 
	"_condition-sex_rev-gene-input_covars-pc3-age-nspots.rda"))
#head(coef(results))

comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")

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

ggsave(paste0("plots/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set, 
	"_rev-gene-input_p-val-histogram-volcano.png"),
	gridExtra::grid.arrange(sex.res$phist, p, ncol=1),
	bg="white", width=6, height=11)
cat("\n\nSaved un-adjusted p value histogram and volcano plots to:",
	paste0("plots/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set, 
	"_rev-gene-input_p-val-histogram-volcano.png"), "\n")

write.csv(eb.f, paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,
	"_rev-gene-input_F-test.csv"), row.names=F)
cat("\n\n\nSaved F test results dframe to:", paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,
        "_rev-gene-input_F-test.csv"), "\n")
write.csv(sex.res$results, paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set, 
	"_rev-gene-input_moderated-t-test.csv"), row.names=F)
cat("\nSaved moderated t test results dframe to:", paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set, 
	"_rev-gene-input_moderated-t-test.csv"), "\n\n")



### layer restricted
results <- readRDS(paste0("processed-data/07_dx_DE/lmFit-voom_layer-restricted_", results_set, 
	"_condition-sex_rev-gene-input_covars-pc3-age-nspots.rda"))


#need to have "dot" instead of "." in comparisons for cluster name 
comparisons = unlist(lapply(comp_names, 
                            function(x) paste(x, comparisons, sep="_")))

eb <- eBayes(groupContrasts(results, comparisons), trend=T, robust=T)
#plotSA(eb)

## F statistics
eb_f = getTopTable(eb, .coef="all") %>%
	mutate(cluster="all", F_design="~ 0 + dx*sex*cluster")

### nested approach
eb_f.list <- lapply(names(comp_names), function(x) {
  tmp = getTopTable(eb, .coef=grep(x, colnames(coef(eb)), value=T))
  colnames(tmp)[1:6] = sapply(strsplit(colnames(tmp)[1:6], "_"), function(y) paste(y[[2]], y[[3]], sep="_"))
  tmp$cluster = x
  tmp$F_design = paste0(x, ": ~ 0 + dx*sex")
  return(tmp)
  })
### implement global padj method (without this, different t cutoffs were sig because of different group sizes)
### details: https://github.com/cran/limma/blob/master/R/decidetests.R
eb_f.list_df = do.call(rbind, eb_f.list) %>% 
  mutate(cluster=factor(cluster, levels=names(comp_names)))
eb_f.list_df$adj.P.Val_global = p.adjust(eb_f.list_df$P.Value, method="BH")

#f_sig_list = lapply(names(comp_names), function(x)
#	filter(eb_f.list_df, cluster==x, adj.P.Val<.05)$gene_id
#)
#names(f_sig_list) <- names(comp_names)

## moderated t statistics
sex.resList <- lapply(names(comp_names), function(x) sexTopTable(eb, phist=F, cluster=x))
sex.res_df = do.call(rbind, sex.resList) %>%
  mutate(cluster=factor(cluster, levels=names(comp_names)))
### implement global padj method (without this, different t cutoffs were sig because of different group sizes)
### details: https://github.com/cran/limma/blob/master/R/decidetests.R
sex.res_df$adj.P.Val = p.adjust(sex.res_df$P.Value, method="BH")

#alt pvalue histogram
tmp = mutate(sex.res_df, sex.group=factor(paste(group, sex), levels=c("NTC.MDD F","NTC.MDD M","NTC.BPD F","NTC.BPD M","MDD.BPD F","MDD.BPD M")))
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

p1 <- ggplot(filter(sex.res_df, sex=="F"), aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point(size=.3)+geom_point(data=filter(tmp, sex=="F"), size=.3, color="red2")+
  facet_grid(rows=vars(cluster), cols=vars(group))+
  xlim(-ceiling(max(sex.res_df$logFC)), ceiling(max(sex.res_df$logFC)))+ggtitle("Females")+
  theme_bw()

p2 <- ggplot(filter(sex.res_df, sex=="M"), aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point(size=.3)+geom_point(data=filter(tmp, sex=="M"), size=.3, color="red2")+
  facet_grid(rows=vars(cluster), cols=vars(group))+
  xlim(-ceiling(max(sex.res_df$logFC)), ceiling(max(sex.res_df$logFC)))+ggtitle("Males")+
  theme_bw()

pdf(file=paste0("plots/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set, 
	"_rev-gene-input_p-val-histogram-volcano.pdf"),
    width=7, height=8)
p0
ggrastr::rasterize(p1,layer='point',dpi=300)
ggrastr::rasterize(p2,layer='point',dpi=300)
dev.off()
cat("\n\nSaved volcano plots to:", paste0("plots/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set, 
	"_rev-gene-input_volcano.pdf"),"\n\n\n")


ph1 <- ggplot(eb.f, aes(x=P.Value))+
  geom_histogram(bins=30)+scale_x_continuous(breaks=c(0,.5,1))+
  labs(subtitle=paste0("Layer-adjusted (", results_set, ")"),
    title="~ 0 + dx*sex", x="F test p value")+
  theme_minimal()+theme(panel.grid.minor=element_blank())

ph2 <- ggplot(eb_f, aes(x=P.Value))+
  geom_histogram(bins=30)+scale_x_continuous(breaks=c(0,.5,1))+
  labs(subtitle=paste0("Layer-restricted (", results_set, ")"),
    title="~ 0 + dx*sex*cluster", x="F test p value")+
  theme_minimal()+theme(panel.grid.minor=element_blank())

ph3 <- ggplot(eb_f.list_df, aes(x=P.Value))+
  geom_histogram(bins=30)+scale_x_continuous(breaks=c(0,.5,1))+
  facet_wrap(vars(cluster), ncol=2, scales="free_y")+
  labs(subtitle=paste0("Layer-restricted (", results_set, ")"), 
    title="cluster: ~ 0 + dx*sex", x="F test p value")+
  theme_minimal()+theme(panel.grid.minor=element_blank())

lay_mat= rbind(c(1,2),c(3,3),c(3,3))
ggsave(paste0("plots/07_dx_DE/F-test_la-lr-pc3-age-nspots_", results_set,
        "_rev-gene-input_p-val-histogram.png"),
        gridExtra::grid.arrange(ph1, ph2, ph3, layout_matrix=lay_mat),
        bg="white", width=6, height=11)
cat("\n\nSaved F test un-adjusted p value histogram to:",
        paste0("plots/07_dx_DE/F-test_la-lr-pc3-age-nspots_", results_set,
        "_rev-gene-input_p-val-histogram-volcano.png"), "\n")


write.csv(eb_f, paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
	"_rev-gene-input_F-test_all.csv"), row.names=F)
write.csv(eb_f.list_df, paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
        "_rev-gene-input_F-test_each.csv"), row.names=F)
cat("\n\n\nSaved F test results dframe to:", paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
        "_rev-gene-input_F-test.csv"), "\n")
write.csv(sex.res_df, paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
	"_rev-gene-input_moderated-t-test.csv"), row.names=F)
cat("\nSaved moderated t test results dframe to:", paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
	"_rev-gene-input_moderated-t-test.csv"), "\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
