setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
})
set.seed(123)

comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")
names(comparisons) <- comparisons
comparisons2 = comparisons[c(1,3,5,2,4,6)]

covars = c("none","pc3 only","nspots","chrM_ratio","age","BMI","Smoking","RIN")
names(covars) <- covars

# layer adjusted
la.all = do.call(rbind, lapply(c("smoothed-k9-1663","seurat-pc30"), function(results_set) {
  ## F test
  f1 = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-no-covars_", 
                       results_set, "_rev-gene-input_F-test.csv"), row.names=1) %>%
    mutate(covar="none")
  f2 = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3_", 
                       results_set, "_rev-gene-input_F-test.csv")) %>%
    mutate(covar="pc3 only")
  f3 = read.csv(paste0("processed-data/07-1_covariate_sensitivity/layer-adjusted_", 
                       results_set, "_rev-gene-input_F-test_all-covar-models.csv"))
  
  la.all = bind_rows(f1, f2, f3) %>% 
    filter(covar %in% covars) %>%
    mutate(covar=factor(covar, levels=covars), annotation=results_set)
  return(la.all)
})) %>% mutate(annotation=factor(annotation, levels=c("smoothed-k9-1663", "seurat-pc30"), 
                                 labels=c("smoothed","seurat")))

## t test
lat.all = do.call(rbind, lapply(c("smoothed-k9-1663","seurat-pc30"), function(results_set) {
  t1 = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-no-covars_", 
                       results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
    mutate(covar="none")
  
  t2 = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3_", 
                       results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
    mutate(covar="pc3 only")
  
  t3 = read.csv(paste0("processed-data/07-1_covariate_sensitivity/layer-adjusted_", 
                       results_set, "_rev-gene-input_moderated-t-test_all-covar-models.csv"))
  
  la.t.all = bind_rows(t1, t2, t3) %>% 
    filter(covar %in% covars) %>%
    mutate(coef= factor(coef, levels=comparisons),
           covar=factor(covar, levels=covars),
           annotation= results_set)
  return(la.t.all)
})) %>% mutate(annotation=factor(annotation, levels=c("smoothed-k9-1663", "seurat-pc30"), 
                                labels=c("smoothed","seurat")))

### compare number of DEGs
la.degs = do.call(rbind, lapply(covars, function(y) {
  sig.genes = filter(la.all, covar==y, annotation=="smoothed", adj.P.Val<.05)$gene_id
  f1 = filter(lat.all, covar==y, adj.P.Val<.05, annotation=="smoothed", gene_id %in% sig.genes)
  sig.genes = filter(la.all, covar==y, annotation=="seurat", adj.P.Val<.05)$gene_id
  f2 = filter(lat.all, covar==y, adj.P.Val<.05, annotation=="seurat", gene_id %in% sig.genes)
  return(rbind(f1, f2))
})) 


# combine all results for L-A model
df1 = bind_rows(filter(la.all, adj.P.Val<.05) %>% group_by(annotation, covar, .drop=F) %>% tally() %>%
  mutate(statistic="F-test"),
  filter(lat.all, adj.P.Val<.05) %>% group_by(covar, annotation, .drop=F) %>% tally() %>%
    mutate(statistic="t-test"),
  group_by(la.degs, annotation, covar, .drop=F) %>% tally() %>% mutate(statistic="F-test and t-test"))


# layer-restricted
lr.all = do.call(rbind, lapply(c("smoothed-k9-1663","seurat-pc30"), function(results_set) {
  ## F test
  f1.all = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-no-covars_", 
                           results_set, "_rev-gene-input_F-test_all.csv")) %>%
    mutate(covar="none")
  
  f2.all = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3_", 
                           results_set, "_rev-gene-input_F-test_all.csv")) %>%
    mutate(covar="pc3 only")
  
  f3.all = read.csv(paste0("processed-data/07-1_covariate_sensitivity/layer-restricted_", 
                           results_set, "_rev-gene-input_F-test_all_all-covar-models.csv"))
  
  lr.all = bind_rows(f1.all, f2.all, f3.all) %>%
    filter(covar %in% covars) %>%
    mutate(covar=factor(covar, levels=covars), annotation= results_set)
  return(lr.all[,c("AveExpr","F","P.Value","adj.P.Val","df.prior","df.total","gene_id","gene_name",
                   "cluster","F_design","covar","annotation")])
})) %>% mutate(annotation=factor(annotation, levels=c("smoothed-k9-1663", "seurat-pc30"), 
                                labels=c("smoothed","seurat")))

## t test
lrt.all = do.call(rbind, lapply(c("smoothed-k9-1663","seurat-pc30"), function(results_set) {
  if(results_set == "smoothed-k9-1663") clus_levels=c("L1","L2","L3.4","L5","L6","WM")
  if(results_set == "seurat-pc30") clus_levels=c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo")
  
  t1 = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-no-covars_", 
                       results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
    mutate(covar="none")
  
  t2 = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3_", 
                       results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
    mutate(covar="pc3 only")
  
  t3 = read.csv(paste0("processed-data/07-1_covariate_sensitivity/layer-restricted_", 
                       results_set, "_rev-gene-input_moderated-t-test_all-covar-models.csv"))
  
  lr.t.all = bind_rows(t1, t2, t3) %>% 
    filter(covar %in% covars) %>%
    mutate(sex.group= factor(paste(sex, group, sep="_"), levels=comparisons),
           covar=factor(covar, levels=covars),
           cluster=factor(cluster, levels=clus_levels),
           annotation= results_set)
  return(lr.t.all)
  })) %>% mutate(annotation=factor(annotation, levels=c("smoothed-k9-1663", "seurat-pc30"), 
                                   labels=c("smoothed","seurat")))


### compare number of DEGs
lr.degs = do.call(rbind, lapply(covars, function(y) {
  sig.genes = filter(lr.all, covar==y, annotation=="smoothed", adj.P.Val<.05)$gene_id
  f1 = filter(lrt.all, covar==y, adj.P.Val<.05, annotation=="smoothed", gene_id %in% sig.genes)
  sig.genes = filter(lr.all, covar==y, annotation=="seurat", adj.P.Val<.05)$gene_id
  f2 = filter(lrt.all, covar==y, adj.P.Val<.05, annotation=="seurat", gene_id %in% sig.genes)
  return(rbind(f1, f2))
})) 


# combine all results for L-R
df2 = bind_rows(filter(lr.all, adj.P.Val<.05) %>% group_by(annotation, covar, .drop=F) %>% tally() %>%
                  mutate(statistic="F-test"),
                filter(lrt.all, adj.P.Val<.05) %>% group_by(covar, annotation, .drop=F) %>% tally() %>%
                  mutate(statistic="t-test"),
                group_by(lr.degs, annotation, covar, .drop=F) %>% tally() %>% mutate(statistic="F-test and t-test"))


# combine L-A and L-R results
df.both = bind_rows(mutate(df1, model="L-A"), mutate(df2, model="L-R")) %>%
  mutate(model=factor(model, levels=c("L-A","L-R")),
         statistic=factor(statistic, levels=c("F-test","t-test","F-test and t-test")))

p1 <- ggplot(df.both, aes(x=covar, y=n, fill=annotation))+
  geom_bar(stat="identity", position="dodge", color="black", linewidth=.5)+
  scale_fill_manual(values=c("white","grey"))+
  facet_grid(cols=vars(statistic), rows=vars(model))+
  labs(title="gene repeats counted", fill="")+
  theme_minimal()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
	legend.text=element_text(size=6), legend.key.size=unit(6,"pt"), panel.grid.minor=element_blank())


df1_single = bind_rows(filter(la.all, adj.P.Val<.05) %>% group_by(annotation, covar, .drop=F) %>% tally() %>%
                  mutate(statistic="F-test"),
                filter(lat.all, adj.P.Val<.05) %>% distinct(covar, annotation, gene_id) %>% 
                  group_by(covar, annotation, .drop=F) %>% tally() %>%
                  mutate(statistic="t-test"),
                distinct(la.degs, covar, annotation, gene_id) %>% 
                  group_by(annotation, covar, .drop=F) %>% tally() %>% 
                  mutate(statistic="F-test and t-test"))


df2_single = bind_rows(filter(lr.all, adj.P.Val<.05) %>% group_by(annotation, covar, .drop=F) %>% tally() %>%
                  mutate(statistic="F-test"),
                filter(lrt.all, adj.P.Val<.05) %>% distinct(covar, annotation, sex.group, gene_id) %>% 
                  group_by(covar, annotation, .drop=F) %>% tally() %>%
                  mutate(statistic="t-test"),
                distinct(lr.degs, covar, annotation, sex.group, gene_id) %>% 
                  group_by(annotation, covar, .drop=F) %>% tally() %>% 
                  mutate(statistic="F-test and t-test"))

df.both_single = bind_rows(mutate(df1_single, model="L-A"), mutate(df2_single, model="L-R")) %>%
  mutate(model=factor(model, levels=c("L-A","L-R")),
         statistic=factor(statistic, levels=c("F-test","t-test","F-test and t-test")))

p2 <- ggplot(df.both_single, aes(x=covar, y=n, fill=annotation))+
  geom_bar(stat="identity", position="dodge", color="black", linewidth=.5)+
  scale_fill_manual(values=c("white","grey"))+
  facet_grid(cols=vars(statistic), rows=vars(model))+
  labs(title="limit 1 count per gene per dx*sex group", fill="")+
  theme_minimal()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
	legend.text=element_text(size=6), legend.key.size=unit(6,"pt"), panel.grid.minor=element_blank())

pdf(file="plots/publication/supp_covariate-selection/DE-results-per-covar_barplot.pdf", height=4, width=6.5)
p1
p2
dev.off()

# look by dx*sex group
# single gene (no repeats if multiple contrasts)
test1 = distinct(la.degs, covar, annotation, coef, gene_id) %>% 
  group_by(annotation, covar, coef, .drop=F) %>% tally() %>% 
  mutate(coef=factor(coef, levels=comparisons2))

col.pal = c("grey50", RColorBrewer::brewer.pal(n=length(covars)-1, "Dark2"))
names(col.pal) = covars

p1 <- ggplot(test1, aes(x=coef, y=n, color=covar))+
  #ggbeeswarm::geom_beeswarm(aes(shape=annotation), size=2)+
  geom_point(aes(shape=annotation), size=2)+
  scale_shape_manual(values=c(5,1), guide="none")+scale_color_manual(values=col.pal)+
  geom_line(data=filter(test1, annotation=="smoothed"), aes(group=covar), lty=1)+
  geom_line(data=filter(test1, annotation=="seurat"), aes(group=covar), lty=2)+
  scale_x_discrete(labels=gsub("_", "\n", comparisons))+
  labs(title="Layer-adjusted DEGs (t test <.05 and F test <.05)", 
       y="# unique genes", x="")+
  theme_minimal()+theme(legend.text=element_text(size=6), legend.key.size=unit(6,"pt"), panel.grid.minor=element_blank(),
	plot.title=element_text(size=8))

test2 = distinct(lr.degs, covar, annotation, sex.group, gene_id) %>% 
  group_by(annotation, covar, sex.group, .drop=F) %>% tally() %>% 
  mutate(coef=factor(sex.group, levels=comparisons2))

p2 <- ggplot(test2, aes(x=coef, y=n, color=covar))+
  #ggbeeswarm::geom_beeswarm(aes(shape=annotation), size=2)+
  geom_point(aes(shape=annotation), size=2)+
  scale_shape_manual(values=c(5,1), guide="none")+scale_color_manual(values=col.pal)+
  geom_line(data=filter(test2, annotation=="smoothed"), aes(group=covar), lty=1)+
  geom_line(data=filter(test2, annotation=="seurat"), aes(group=covar), lty=2)+
  scale_x_discrete(labels=gsub("_", "\n", comparisons))+
  labs(title="Layer-restricted DEGs (t test <.05 and F test <.05)",
       y="# unique genes", x="")+
  theme_minimal()+theme(legend.text=element_text(size=6), legend.key.size=unit(6,"pt"), panel.grid.minor=element_blank(),
	plot.title=element_text(size=8))

# add final model results
source("code/09_DEG_GRN/load_DEGs.r")
la.sig.df = bind_rows(group_by(sigList[["sm.la"]], sex.group, .drop=F) %>% tally() %>%
                        mutate(annotation="smoothed", model="L-A"),
                      group_by(sigList[["se.la"]], sex.group, .drop=F) %>% tally() %>%
                        mutate(annotation="seurat", model="L-A")) %>%
  mutate(annotation=factor(annotation, levels=c("smoothed","seurat")),
         model=factor(model, levels=c("L-A","L-R")),
         covar="final model")
lr.sig.df = bind_rows(distinct(sigList[["sm.lr"]], sex.group, gene_id) %>% 
                        group_by(sex.group, .drop=F) %>% tally() %>%
                        mutate(annotation="smoothed", model="L-R"),
                      distinct(sigList[["se.lr"]], sex.group, gene_id) %>%
                        group_by(sex.group, .drop=F) %>% tally() %>%
                        mutate(annotation="seurat", model="L-R")) %>%
  mutate(annotation=factor(annotation, levels=c("smoothed","seurat")),
         model=factor(model, levels=c("L-A","L-R")),
         covar="final model")

p1_mod <- p1+geom_line(data=filter(la.sig.df, annotation=="smoothed"), aes(x=sex.group, group=covar), linewidth=1, lty=1, color="black")+
  geom_line(data=filter(la.sig.df, annotation=="seurat"), aes(x=sex.group, group=covar), linewidth=1, lty=2, color="black")
p2_mod <- p2+geom_line(data=filter(lr.sig.df, annotation=="smoothed"), aes(x=sex.group, group=covar), linewidth=1, lty=1, color="black")+
  geom_line(data=filter(lr.sig.df, annotation=="seurat"), aes(x=sex.group, group=covar), linewidth=1, lty=2, color="black")

pdf(file="plots/publication/supp_covariate-selection/DE-results-per-covar-dx-sex_lineplot.pdf", height=4, width=6.5)
grid.arrange(p1_mod+ylim(0,200), p2_mod+ylim(0,200), ncol=1)
dev.off()


cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
