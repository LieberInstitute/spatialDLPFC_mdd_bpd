setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
})

DATASET_ID = "spe-n119_13844-no-lowUMI_adj_with-logcounts-corr"

# load in GRNBOOST2 adj output
lg.mask = read.csv(paste0("processed-data/09_DEG_GRN/", DATASET_ID, ".csv"))


# load avg expr
avg.expr = read.csv("processed-data/06_pseudobulk/PRECAST_smoothed/pseudobulk-sample-smoothed-n1663-k9_filtered-genes_avg-logcounts.csv", row.names=1)
#avg.expr = read.csv("processed-data/06_pseudobulk/Seurat/pseudobulk-sample-seurat-pc30-no-lowUMI_filtered-genes_avg-logcounts.csv", row.names=1)

# load rowData input to filter avg expr to input for GRN
load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
#load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-no-lowUMI_sample-seurat-pc30_norm-filt.Rdata")

rowData(spe_pseudo)$high_expr_group_sample_id2 <- edgeR::filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster2 <- edgeR::filterByExpr(spe_pseudo, group = spe_pseudo$smoothed_k9_1663)
#rowData(spe_pseudo)$high_expr_group_cluster2 <- edgeR::filterByExpr(spe_pseudo, group = spe_pseudo$seurat_label)

rdata = as.data.frame(rowData(spe_pseudo)[,c("gene_id","gene_name","gene_type")])
rdata$DE_input = rowData(spe_pseudo)$high_expr_group_sample_id2==T & rowData(spe_pseudo)$high_expr_group_cluster2==T

# rename gene with duplicated gene name
avg.expr = avg.expr[rdata$gene_id[rdata$DE_input],]

# precast smoothed
avg.expr[c("ENSG00000187522","ENSG00000271858"),"gene_name"] = c("MSTANTD7","LOC101928965")

## seurat labels
##table(rdata$DE_input) #only 1 gene
#avg.expr[c("ENSG00000187522"),"gene_name"] = "MSTANTD7"

decile.pal = rainbow(10)
names(decile.pal) = paste0("dec",1:10)

# filters to top 19 because creating modules automatically add TF self to module
cat("\nModule criteria: abs(rho)>.03, importance>.5, >=20 target genes\n")
modules = filter(lg.mask, importance>.5, regulation!=0) %>% 
  group_by(regulation, TF) %>% add_tally() %>%
  filter(n>=19) 
#modules$n <- NULL

act_mod = nrow(filter(distinct(modules, TF, regulation, n), regulation>0))
rep_mod = nrow(filter(distinct(modules, TF, regulation, n), regulation<0))
cat("\nNumber of activating modules matching criteria:", act_mod,"\n")
cat("\nNumber of activating modules matching criteria:", rep_mod, "\n")

cat("\nModule criteria: abs(rho)>.03, top 20% most importance target genes per TF, importance>.5, >=20 target genes\n")
modules2 = filter(lg.mask, regulation!=0) %>% 
  group_by(regulation, TF) %>% slice_max(prop = .2, importance) %>%
  filter(importance>.5) %>% 
  group_by(TF, regulation) %>% add_tally() %>%
  filter(n>=19)
act_mod2 = nrow(filter(distinct(modules2, TF, regulation, n), regulation>0))
rep_mod2 = nrow(filter(distinct(modules2, TF, regulation, n), regulation<0))
cat("\nNumber of activating modules matching criteria:", act_mod2,"\n")
cat("\nNumber of activating modules matching criteria:", rep_mod2, "\n")


t1 = group_by(modules, TF, dir=paste0("rho_",factor(regulation, levels=c(-1,1), labels=c("rep","act")))) %>% 
  tally()
t2 = group_by(modules2, TF, dir=paste0("rho_",factor(regulation, levels=c(-1,1), labels=c("rep","act")))) %>% 
  tally()

## use to set axis limits
cat("\nMax number of genes for DEG modules:", max(c(t1$n, t2$n)), "\n")

t3 = full_join(t1, t2, by=c("TF","dir"), suffix=c("_not20","_top20")) %>%
  left_join(select(avg.expr, gene_name, TF_decile=decile), by=c("TF"="gene_name")) %>%
  mutate(TF_decile=factor(TF_decile, levels=c(1:10), labels=paste0("dec",1:10)),
         dir=factor(dir, levels=c("rho_act","rho_rep"), labels=c("Activating module","Repressing module")))

plist = list()
plist[[1]] = ggplot(t3, aes(x=n_not20, y=n_top20, color=TF_decile))+
  geom_point()+
  scale_x_log10("# target genes", limits=c(10,13000))+
  scale_y_log10("# target genes (top 20% filter)", limits=c(10,13000))+
  scale_color_manual("TF expr.\ndecile", values=decile.pal)+
  facet_wrap(vars(dir), ncol=1)+
  labs(title="Number of size of TF-target gene modules", subtitle=DATASET_ID)



# look at dominance of highly expressed genes in modules

t1 = left_join(modules, select(avg.expr, gene_name, "decile_target"=decile), by=c("target"="gene_name")) %>%
  mutate(decile_target=factor(decile_target, levels=1:10)) %>%
  group_by(TF, regulation, n, decile_target, .drop=F) %>% tally(name="n_decile") %>%
  mutate(prop_decile=n_decile/n)

p1 <- ggplot(t1, aes(x=decile_target, y=prop_decile))+
  ggbeeswarm::geom_quasirandom(width = .2)+
  geom_boxplot(alpha=.5, color="red3", width=.5, outliers=F)+
  labs(title="Prop. module targets in expr. deciles",
       y="prop. of module target genes", x="target gene expr. decile", subtitle=DATASET_ID)

t2 = left_join(modules2, select(avg.expr, gene_name, "decile_target"=decile), by=c("target"="gene_name")) %>%
  mutate(decile_target=factor(decile_target, levels=1:10)) %>%
  group_by(TF, regulation, n, decile_target, .drop=F) %>% tally(name="n_decile") %>%
  mutate(prop_decile=n_decile/n)

p2 <- ggplot(t2, aes(x=decile_target, y=prop_decile))+
  ggbeeswarm::geom_quasirandom(width = .2)+
  geom_boxplot(alpha=.5, color="red3", width=.5, outliers=F)+
  labs(title="Prop. module targets in expr. deciles (top 20% filter)",
       y="prop. of module target genes", x="target gene expr. decile", subtitle=DATASET_ID)

plist[[2]] <- grid.arrange(p1, p2, ncol=1)
  


t1.1 = left_join(t1, select(avg.expr, gene_name, "decile_TF"=decile), by=c("TF"="gene_name")) %>%
  mutate(decile_TF=factor(decile_TF, levels=10:1, labels=paste0("TF_dec",10:1)))
t1.1$decile_TF = droplevels(t1.1$decile_TF)

plist[[3]] <- ggplot(t1.1, aes(x=decile_target, y=prop_decile))+
  ggbeeswarm::geom_quasirandom(size=.5, width = .2)+
  geom_boxplot(alpha=.5, color="red3", width=.5, outliers=F)+
  facet_wrap(vars(decile_TF), ncol=1)+
  labs(title="Prop. module targets in expr. deciles (facet by TF decile)",
       y="prop. of module target genes", x="target gene expr. decile", subtitle=DATASET_ID)


t2.1 = left_join(t2, select(avg.expr, gene_name, "decile_TF"=decile), by=c("TF"="gene_name")) %>%
  mutate(decile_TF=factor(decile_TF, levels=10:1, labels=paste0("TF_dec",10:1)))
t2.1$decile_TF = droplevels(t2.1$decile_TF)

plist[[4]] <- ggplot(t2.1, aes(x=decile_target, y=prop_decile))+
  ggbeeswarm::geom_quasirandom(size=.5, width = .2)+
  geom_boxplot(alpha=.5, color="red3", width=.5, outliers=F)+
  facet_wrap(vars(decile_TF), ncol=1)+
  labs(title="Prop. module targets in expr. deciles (top 20% filter) (facet by TF decile)",
       y="prop. of module target genes", x="target gene expr. decile", subtitle=DATASET_ID)



ggsave(file=paste0("plots/09_DEG_GRN/", DATASET_ID, "_explore-modules.pdf"), 
       marrangeGrob(grobs=plist, ncol=1, nrow=1, top=NULL))
cat("\nPlots saved to:", paste0("plots/09_DEG_GRN/", DATASET_ID, "_explore-modules.pdf"), "\n")

modules2$n <- NULL
write.csv(modules2, paste0("processed-data/09_DEG_GRN/", DATASET_ID, "_top20percent.csv"), row.names=F)
cat("\nFiltered adj. list for top 20% modules saved to:", paste0("processed-data/09_DEG_GRN/", DATASET_ID, "_top20percent.csv"), "\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
