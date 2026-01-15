library(dplyr)
library(ggplot2)

results_set = "smoothed-k9-1663"
comp_names = c("L1","L2","L3dot4","L5","L6","WM")
names(comp_names) = c("L1","L2","L3.4","L5","L6","WM")


results_set = "seurat-pc30"
comp_names = c("MicrodotVasc","Astro","L2dot3","L4","Inhb","L5","L6","Oligo")
names(comp_names) = c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo")


comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")
names(comparisons) <- comparisons

#model with pc3 only ----

covars = c("none","pc3 only","age","BMI","Smoking","RIN","PMI")
names(covars) <- covars

## layer adjusted ----
f1 = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-no-covars_", 
                     results_set, "_rev-gene-input_F-test.csv"), row.names=1) %>%
  mutate(covar="none")
f2 = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3_", 
                     results_set, "_rev-gene-input_F-test.csv")) %>%
  mutate(covar="pc3 only")
f3 = read.csv(paste0("processed-data/07-1_covariate_sensitivity/layer-adjusted_", 
                     results_set, "_rev-gene-input_F-test_all-covar-models.csv"))

la.all = bind_rows(f1, f2, f3) %>% 
  mutate(covar=factor(covar, levels=covars))

filter(la.all, adj.P.Val<.05) %>% group_by(covar) %>% tally()
filter(la.all, adj.P.Val<.05) %>% group_by(covar) %>% slice_min(F) %>%
  select(covar, F, P.Value, adj.P.Val)

la.all.baseline = filter(la.all, adj.P.Val<.05, covar=="none") %>% group_by(covar) %>% tally(name="n_none") 
la.all2 = filter(la.all, adj.P.Val<.05) %>% group_by(covar) %>% tally(name="n_covar_model") %>%
  mutate(n_none=la.all.baseline$n_none) %>% filter(covar!="none") %>% mutate(prop_no.covar_degs=n_covar_model/n_none)


## layer-restricted ----
f1 = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-no-covars_", 
                     results_set, "_rev-gene-input_F-test.csv"))
#smoothed
f1.all = mutate(f1[is.na(f1$F_NTC.MDD),1:46], covar="none")
f1.each = mutate(f1[!is.na(f1$F_NTC.MDD),c(47:52,37:46)], covar="none")
#seurat
f1.all = mutate(f1[is.na(f1$F_NTC.MDD),1:58], covar="none")
f1.each = mutate(f1[!is.na(f1$F_NTC.MDD),c(59:64,49:58)], covar="none")


f2 = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3_", 
                     results_set, "_rev-gene-input_F-test.csv"))
#smoothed
f2.all = mutate(f2[is.na(f2$F_NTC.MDD),1:46], covar="pc3 only")
f2.each = mutate(f2[!is.na(f2$F_NTC.MDD),c(47:52,37:46)], covar="pc3 only")
#seurat
f2.all = mutate(f2[is.na(f2$F_NTC.MDD),1:58], covar="pc3 only")
f2.each = mutate(f2[!is.na(f2$F_NTC.MDD),c(59:64,49:58)], covar="pc3 only")

f3 = read.csv(paste0("processed-data/07-1_covariate_sensitivity/layer-restricted_", 
                     results_set, "_rev-gene-input_F-test_all-covar-models.csv"))
#smoothed
f3.all = f3[is.na(f3$F_NTC.MDD),1:47]
f3.each = f3[!is.na(f3$F_NTC.MDD),c(48:53,37:47)]
#seurat
f3.all = f3[is.na(f3$F_NTC.MDD),1:59]
f3.each = f3[!is.na(f3$F_NTC.MDD),c(60:65,49:59)]

lr.all = bind_rows(f1.all, f2.all, f3.all) %>%
  mutate(covar=factor(covar, levels=covars))


filter(lr.all, adj.P.Val<.05) %>% group_by(covar) %>% tally()
filter(lr.all, adj.P.Val<.05) %>% group_by(covar) %>% slice_min(F) %>%
  select(covar, F, P.Value, adj.P.Val)

lr.all.baseline = filter(lr.all, adj.P.Val<.05, covar=="none") %>% group_by(cluster, .drop=F) %>% tally(name="n_none") 
lr.all2 = filter(lr.all, adj.P.Val<.05) %>% group_by(covar) %>% tally(name="n_covar_model") %>%
  mutate(n_none=lr.all.baseline$n_none) %>% filter(covar!="none") %>% mutate(prop_no.covar_degs=n_covar_model/n_none)


checkList <- lapply(c("pc3 only", "age","BMI","PMI"), function(x) filter(lr.all, adj.P.Val<.05, covar==x)$gene_id)
names(checkList) <- c("pc3","age","BMI","PMI")

UpSetR::upset(UpSetR::fromList(checkList), text.scale=2)
#lr.all$adj.P.Val_global = p.adjust(lr.all$P.Value, method="BH")
#filter(lr.all, adj.P.Val_global<.05) %>% group_by(covar) %>% tally()
#filter(lr.all, adj.P.Val_global<.05) %>% group_by(covar) %>% slice_min(F) %>%
#  select(covar, F, P.Value, adj.P.Val, adj.P.Val_global)


lr.each = bind_rows(f1.each, f2.each, f3.each) %>%
  mutate(covar=factor(covar, levels=covars),
         cluster=factor(cluster, levels=names(comp_names)))

filter(lr.each, adj.P.Val<.05) %>% group_by(covar, cluster, .drop=F) %>% tally() %>%
  tidyr::pivot_wider(names_from="covar", values_from="n")

filter(lr.each, adj.P.Val<.05) %>% group_by(covar, cluster) %>% slice_min(F) %>%
  select(covar, cluster, F, P.Value, adj.P.Val)


lr.each.baseline = filter(lr.each, adj.P.Val<.05, covar=="none") %>% group_by(cluster, .drop=F) %>% tally(name="n_none") 
lr.each2 = filter(lr.each, adj.P.Val<.05, covar!="none") %>% group_by(covar, cluster) %>% tally(name="n_covar_model")
lr.each2 = left_join(lr.each2, lr.each.baseline) %>% mutate(prop_no.covar_degs=n_covar_model/n_none)

lr.any = bind_rows(lr.each2, mutate(lr.all2, cluster="all"),
                   mutate(la.all2, cluster="L-A")) %>% mutate(cluster= factor(cluster, levels=c("L-A","all",names(comp_names))))

covars.colors = c("white","black","#E4775D","#A0C255","#EEBC4A","grey","#93D3F6")
names(covars.colors) = names(covars)
ggplot(lr.any, aes(x=cluster, y=prop_no.covar_degs, color=covar))+
  ggbeeswarm::geom_beeswarm(size=2)+scale_color_manual(values=covars.colors)+
  geom_hline(aes(yintercept=1), lty=2)+
  labs(title="F test sig. genes")





checkList <- lapply(names(comp_names), function(y) {
  tmp = lapply(c("pc3 only","age","BMI","PMI"), function(x) filter(lr.each, adj.P.Val<.05, cluster==y, covar==x)$gene_id)
  names(tmp) <- c("pc3","age","BMI","PMI")
  return(tmp)
})
names(checkList) <- names(comp_names)
UpSetR::upset(UpSetR::fromList(checkList[["L1"]]), text.scale=2, sets.x.label="L1 Set Size")
UpSetR::upset(UpSetR::fromList(checkList[["L2"]]), text.scale=2, sets.x.label="L2 Set Size")
UpSetR::upset(UpSetR::fromList(checkList[["L3.4"]]), text.scale=2, sets.x.label="L3.4 Set Size")
UpSetR::upset(UpSetR::fromList(checkList[["L5"]]), text.scale=2, sets.x.label="L5 Set Size")
UpSetR::upset(UpSetR::fromList(checkList[["L6"]]), text.scale=2, sets.x.label="L6 Set Size")
UpSetR::upset(UpSetR::fromList(checkList[["WM"]]), text.scale=2, sets.x.label="WM Set Size")


checkList2  <- lapply(c("pc3 only","age","BMI","PMI"), function(y) {
  tmp = lapply(names(comp_names), function(x) filter(lr.each, adj.P.Val<.05, cluster==x, covar==y)$gene_id)
  names(tmp) <- names(comp_names)
  return(tmp)
})
names(checkList2) <- c("pc3","age","BMI","PMI")

UpSetR::upset(UpSetR::fromList(checkList2[["pc3"]]), text.scale=2, nsets=6, sets.x.label="pc3 Set Size")
UpSetR::upset(UpSetR::fromList(checkList2[["age"]]), text.scale=2, nsets=6, sets.x.label="age Set Size")

t3 = read.csv("processed-data/07-1_covariate_sensitivity/layer-restricted_smoothed-k9-1663_rev-gene-input_moderated-t-test_all-covar-models.csv")
head(t3)

new.lr_age = filter(t3, covar=="age") %>% mutate(sex.group=factor(paste(sex, group, sep="_"), levels=comparisons),
                                                 cluster=factor(cluster, levels=names(comp_names)))
do.call(rbind, lapply(names(comp_names), function(x) filter(new.lr_age, adj.P.Val<.05, cluster==x, 
                                                            gene_id %in% checkList[[x]]$age) %>% 
         group_by(sex.group, .drop=F) %>% tally() %>% mutate(cluster=x))) %>%
  tidyr::pivot_wider(names_from="cluster", values_from="n")



new_lr.pmi = filter(t3, covar=="PMI") %>% mutate(sex.group=factor(paste(sex, group, sep="_"), levels=comparisons),
                                                 cluster=factor(cluster, levels=names(comp_names)))

do.call(rbind, lapply(names(comp_names), function(x) filter(new_lr.pmi, adj.P.Val<.05, cluster==x, 
                                                            gene_id %in% checkList[[x]]$PMI) %>% 
                        group_by(sex.group, .drop=F) %>% tally() %>% mutate(cluster=x))) %>%
  tidyr::pivot_wider(names_from="cluster", values_from="n")

filter(new_lr.pmi, sex.group=="M_NTC.MDD", adj.P.Val<.05, cluster=="L1")
filter(new_lr.pmi, sex.group=="M_NTC.MDD", adj.P.Val<.05, !cluster %in% c("L1","WM"))

filter(new_lr.pmi, gene_name=="CX3CR1", group=="NTC.BPD", adj.P.Val<.05)


new_la.pmi = read.csv("processed-data/07-1_covariate_sensitivity/layer-adjusted_smoothed-k9-1663_rev-gene-input_moderated-t-test_all-covar-models.csv") %>%
  filter(covar=="PMI")
filter(new_la.pmi, adj.P.Val<.05, gene_name=="SST")
filter(new_la.pmi, adj.P.Val<.05, gene_name=="CORT")
filter(new_la.pmi, adj.P.Val<.05, gene_name=="SURF1")
filter(new_la.pmi, adj.P.Val<.05, gene_name=="FA2H")
filter(new_la.pmi, adj.P.Val<.05, gene_name=="CLDN11")


filter(new_la.pmi, adj.P.Val<.05, gene_name=="CEBPD")
filter(new_la.pmi, adj.P.Val<.05, gene_name=="ANGPTL4")
filter(new_lr.pmi, adj.P.Val<.05, gene_name=="ANGPTL4")
filter(new_la.pmi, adj.P.Val<.05, gene_name=="APOLD1")
filter(new_lr.pmi, adj.P.Val<.05, gene_name=="APOLD1")
