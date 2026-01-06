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

## layer agnostic
results <- readRDS(paste0("processed-data/07_dx_DE/lmFit-voom_layer-agnostic_", results_set,
        "_condition-sex_rev-gene-input_covars-none.rda"))
#head(coef(results))

comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")

eb <- eBayes(groupContrasts(results, comparisons), trend=T)
#plotSA(eb)

sex.res = sexTopTable(eb, phist=T)
#sex.res$phist #awesome, flat not reverse


cat("\nLayer-agnostic analysis with adj p<.05 and abs(logFC)>.2:\n\n")
filter(sex.res$results, adj.P.Val<.05, abs(logFC)>.2) %>% group_by(sex, group, .drop=F) %>% tally()

p <- ggplot(sex.res$results, aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point(size=.3)+geom_point(data=filter(sex.res$results, adj.P.Val<.05, abs(logFC)>.2), size=.3, color="red2")+
  facet_grid(rows=vars(group), cols=vars(sex))+
  xlim(-ceiling(max(sex.res$results$logFC)), ceiling(max(sex.res$results$logFC)))+
  ggtitle(paste0("Layer-agnostic (", results_set, ")"))+
  theme_bw()

ggsave(paste0("plots/07_dx_DE/layer-agnostic-no-covars_", results_set,
        "_rev-gene-input_p-val-histogram-volcano.png"),
        gridExtra::grid.arrange(sex.res$phist, p, ncol=1),
        bg="white", width=6, height=11)
cat("\n\nSaved un-adjusted p value histogram and volcano plots to:",
        paste0("plots/07_dx_DE/layer-agnostic-no-covars_", results_set,
        "_rev-gene-input_p-val-histogram-volcano.png"), "\n")

write.csv(sex.res$results, paste0("processed-data/07_dx_DE/layer-agnostic-no-covars_", results_set,
        "_rev-gene-input_compiled-results.csv"))
cat("\n\n\nSaved compiled results dframe to:", paste0("processed-data/07_dx_DE/layer-agnostic-no-covars_", results_set,
        "_rev-gene-input_compiled-results.csv"), "\n\n")

## layer adjusted
results <- readRDS(paste0("processed-data/07_dx_DE/lmFit-voom_layer-adjusted_", results_set, 
	"_condition-sex_rev-gene-input_covars-none.rda"))
#head(coef(results))

comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")

eb <- eBayes(groupContrasts(results, comparisons), trend=T)
#plotSA(eb)

sex.res = sexTopTable(eb, phist=T)
#sex.res$phist #awesome, flat not reverse

cat("\nLayer-adjusted analysis with adj p<.05 and abs(logFC)>.2:\n\n")
filter(sex.res$results, adj.P.Val<.05, abs(logFC)>.2) %>% group_by(sex, group, .drop=F) %>% tally()

p <- ggplot(sex.res$results, aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point(size=.3)+geom_point(data=filter(sex.res$results, adj.P.Val<.05, abs(logFC)>.2), size=.3, color="red2")+
  facet_grid(rows=vars(group), cols=vars(sex))+
  xlim(-ceiling(max(sex.res$results$logFC)), ceiling(max(sex.res$results$logFC)))+
  ggtitle(paste0("Layer-adjusted (", results_set, ")"))+
  theme_bw()

ggsave(paste0("plots/07_dx_DE/layer-adjusted-no-covars_", results_set, 
	"_rev-gene-input_p-val-histogram-volcano.png"),
	gridExtra::grid.arrange(sex.res$phist, p, ncol=1),
	bg="white", width=6, height=11)
cat("\n\nSaved un-adjusted p value histogram and volcano plots to:",
	paste0("plots/07_dx_DE/layer-adjusted-no-covars_", results_set, 
	"_rev-gene-input_p-val-histogram-volcano.png"), "\n")

write.csv(sex.res$results, paste0("processed-data/07_dx_DE/layer-adjusted-no-covars_", results_set, 
	"_rev-gene-input_compiled-results.csv"))
cat("\n\n\nSaved compiled results dframe to:", paste0("processed-data/07_dx_DE/layer-adjusted-no-covars_", results_set, 
	"_rev-gene-input_compiled-results.csv"), "\n\n")

### layer restricted
results <- readRDS(paste0("processed-data/07_dx_DE/lmFit-voom_layer-restricted_", results_set, 
	"_condition-sex_rev-gene-input_covars-none.rda"))


#need to have "dot" instead of "." in comparisons for cluster name 
comparisons = unlist(lapply(comp_names, 
                            function(x) paste(x, comparisons, sep="_")))

eb <- eBayes(groupContrasts(results, comparisons), trend=T)
#plotSA(eb)

#get top table for all comparisons
sex.resList <- lapply(names(comp_names), function(x) sexTopTable(eb, phist=T, cluster=x))

##plot pval histogram for all comparisons
#pdf(file=paste0("plots/07_dx_DE/layer-restricted-pc3-only_", results_set, 
#	"_rev-gene-input_p-val-histogram.pdf"),
#    width=5, height=6)
#for(i in 1:length(sex.resList)) {
#	plot(sex.resList[[i]]$phist+ggtitle(names(comp_names)[i]))
#}
#dev.off()
#cat("\n\nSaved un-adjusted p value histogram to:", paste0("plots/07_dx_DE/layer-restricted-pc3-only_", results_set, 
#	"_rev-gene-input_p-val-histogram.pdf"),"\n")

sex.res_df = do.call(rbind, lapply(sex.resList, function(x) x$results)) %>%
  mutate(cluster=factor(cluster, levels=names(comp_names)))

#alt pvalue histogram
tmp = mutate(sex.res_df, sex.group=factor(paste(group, sex), levels=c("NTC.MDD F","NTC.MDD M","NTC.BPD F","NTC.BPD M","MDD.BPD F","MDD.BPD M")))
p0 <- ggplot(tmp, aes(x=P.Value))+
  geom_histogram(bins=30)+scale_x_continuous(breaks=c(0,.5,1))+
  facet_grid(rows=vars(cluster), cols=vars(sex.group), scales="free_y")+
  labs(title=paste0("Layer-restricted (", results_set, ")"))+
  theme_minimal()+theme(panel.grid.minor=element_blank())

#volcanos
p1 <- ggplot(filter(sex.res_df, sex=="F"), aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point(size=.3)+geom_point(data=filter(sex.res_df, sex=="F", adj.P.Val<.05, abs(logFC)>.2), size=.3, color="red2")+
  facet_grid(rows=vars(cluster), cols=vars(group))+
  xlim(-ceiling(max(sex.res_df$logFC)), ceiling(max(sex.res_df$logFC)))+ggtitle("Females")+
  theme_bw()

p2 <- ggplot(filter(sex.res_df, sex=="M"), aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point(size=.3)+geom_point(data=filter(sex.res_df, sex=="M", adj.P.Val<.05, abs(logFC)>.2), size=.3, color="red2")+
  facet_grid(rows=vars(cluster), cols=vars(group))+
  xlim(-ceiling(max(sex.res_df$logFC)), ceiling(max(sex.res_df$logFC)))+ggtitle("Males")+
  theme_bw()

pdf(file=paste0("plots/07_dx_DE/layer-restricted-no-covars_", results_set, 
	"_rev-gene-input_p-val-histogram-volcano.pdf"),
    width=7, height=8)
p0
ggrastr::rasterize(p1,layer='point',dpi=300)
ggrastr::rasterize(p2,layer='point',dpi=300)
dev.off()
cat("\n\nSaved volcano plots to:", paste0("plots/07_dx_DE/layer-restricted-no-covars_", results_set, 
	"_rev-gene-input_volcano.pdf"),"\n\n\n")

cat("\nLayer-restricted analysis with adj p<.05 and abs(logFC)>.2:\n\n")
filter(sex.res_df, adj.P.Val<.05, abs(logFC)>.2) %>% group_by(sex, group, cluster, .drop=F) %>% tally() %>%
  tidyr::pivot_wider(names_from="cluster", values_from="n", values_fill=0)

write.csv(sex.res_df, paste0("processed-data/07_dx_DE/layer-restricted-no-covars_", results_set,
	"_rev-gene-input_compiled-results.csv"), row.names=F)
cat("\n\n\nSaved compiled results dframe to:", paste0("processed-data/07_dx_DE/layer-restricted-no-covars_", results_set,
	"_rev-gene-input_compiled-results.csv"), "\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
