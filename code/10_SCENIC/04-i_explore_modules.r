library(dplyr)
library(ggplot2)

results_set="smoothed-k9-1663"
la.degs = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))
lr.degs = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))
mbv.degs = union(la.degs$gene_name, lr.degs$gene_name)
length(mbv.degs) #649

mbv.degs2 = union(la.degs$gene_name[la.degs$n_ttest_sig>0], lr.degs$gene_name[lr.degs$n_ttest_sig>0])
length(mbv.degs2) #387

avg.expr = read.csv("processed-data/06_pseudobulk/PRECAST_smoothed/pseudobulk-sample-smoothed-n1663-k9_filtered-genes_avg-logcounts.csv", row.names=1)
#update genes that had repeated gene names
avg.expr[c("ENSG00000187522","ENSG00000271858","ENSG00000261186","ENSG00000285053","ENSG00000269226"),"gene_name"] = c("MSTANTD7","LOC101928965","LINC03100","GGPS1-TBCE","TMSB15C")
avg.expr = avg.expr[setdiff(rownames(avg.expr), c("ENSG00000286237","ENSG00000261480","ENSG00000280987")),]


# MBv F-test and DEGs by avg.expr decile
ggplot(mutate(avg.expr[union(la.degs$gene_id, lr.degs$gene_id),], 
              is_DEG=factor(as.character(gene_name %in% mbv.degs2), levels=c("FALSE","TRUE"), 
                            labels=c("F-test only","F-test and t-test"))), 
       aes(x=decile, fill=is_DEG))+
  geom_bar(stat="count", position="dodge")+
  scale_fill_manual(values=c("grey70","grey30"))+
  theme_bw()



#logcounts, mask dropouts
lg.mask = read.csv("processed-data/09_SCENIC/spe-n119_13844-no-lowUMI_adj_with-logcounts-corr.csv")

#logcounts, no mask dropouts
lg.nomask = read.csv("processed-data/09_SCENIC/spe-n119_13844-no-lowUMI_adj_with-logcounts-corr-no-mask.csv")


##### look at number of activating and repressing results #####

t1 = group_by(lg.mask, TF, dir=paste0("rho_",factor(sign(rho), levels=c(-1,1), labels=c("rep","act")))) %>% tally() %>%
  tidyr::pivot_wider(names_from="dir", values_from="n", values_fill=0) %>%
  left_join(select(avg.expr, gene_name, "decile_TF"=decile), by=c("TF"="gene_name"))

p1 <- ggplot(t1, aes(x=rho_act, y=rho_rep, color=decile_TF))+
  geom_point(size=.5)+#facet_wrap(vars(decile))+
  scale_color_gradientn(colors=rainbow(10))+
  scale_x_log10("# act. targets (log10 scale)")+
  labs(y="# of rep. targets", title="mask dropouts", subtitle="13844-no-lowUMI_adj_with-logcounts-corr")
       

t2 = group_by(lg.nomask, TF, dir=paste0("rho_",factor(sign(rho), levels=c(-1,1), labels=c("rep","act")))) %>% tally() %>%
  tidyr::pivot_wider(names_from="dir", values_from="n", values_fill=0) %>%
  left_join(select(avg.expr, gene_name, "decile_TF"=decile), by=c("TF"="gene_name"))

p2 <- ggplot(t2, aes(x=rho_act, y=rho_rep, color=decile_TF))+
  geom_point(size=.5)+#facet_wrap(vars(decile))+
  scale_color_gradientn(colors=rainbow(10))+
  scale_x_log10("# act. targets (log10 scale)")+
  labs(y="# of rep. targets", title="no mask dropouts", subtitle="13844-no-lowUMI_adj_with-logcounts-corr")

gridExtra::grid.arrange(p1, p2)

##### check individual genes #####

gene1= "MEF2C" #MEF2C, ELK1
gene2="SOD1" #SOD1, CEBPD
ch1.1 = filter(lg.mask, TF==gene1) 
ch2.1 = filter(lg.nomask, TF==gene1)

ch1.2 = filter(lg.mask, TF==gene2) #CEBPD
ch2.2 = filter(lg.nomask, TF==gene2)


t3.1 = left_join(ch1.1, ch2.1, by=c("TF","target"), suffix=c("_mask","_nomask")) %>%
  mutate(mask_na=is.na(rho_mask)) %>%
  left_join(select(avg.expr, gene_name, "decile_target"=decile), by=c("target"="gene_name"))

p1 <- ggplot(t3.1, aes(x=rho_mask, y=rho_nomask, color=decile_target))+
  geom_point(size=.5)+scale_color_gradientn(colors=rainbow(10))+
  geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+
  ggtitle(gene1)#+labs(x="raw counts", y="logcounts")

t3.2 = left_join(ch1.2, ch2.2, by=c("TF","target"), suffix=c("_mask","_nomask")) %>%
  mutate(mask_na=is.na(rho_mask)) %>%
  left_join(select(avg.expr, gene_name, "decile_target"=decile), by=c("target"="gene_name"))

p2 <- ggplot(t3.2, aes(x=rho_mask, y=rho_nomask, color=decile_target))+
  geom_point(size=.5)+scale_color_gradientn(colors=rainbow(10))+
  geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+
  ggtitle(gene2)#+labs(x="raw counts", y="logcounts")

gridExtra::grid.arrange(p1, p2, ncol=1)


##### THIS IS THE KEY PLOTS #####
# looks at my module filters and classifies the distribution of target genes for each module based on avg expr

t1 = filter(lg.mask, abs(rho)>.03, importance>.5) %>% 
  left_join(select(avg.expr, gene_name, "decile_target"=decile), by=c("target"="gene_name")) %>%
  group_by(TF, regulation) %>% add_tally(name="n_total") %>%
  filter(n_total>=20) %>% mutate(decile_target=as.factor(decile_target)) %>%
  group_by(TF, regulation, n_total, decile_target, .drop=F) %>% tally(name="n_decile") %>%
  mutate(prop_decile=n_decile/n_total)
head(t1)
dim(t1)


t2 = filter(lg.nomask, abs(rho)>.03, importance>.5) %>% 
  left_join(select(avg.expr, gene_name, "decile_target"=decile), by=c("target"="gene_name")) %>%
  group_by(TF, regulation) %>% add_tally(name="n_total") %>%
  filter(n_total>=20) %>% mutate(decile_target=as.factor(decile_target)) %>%
  group_by(TF, regulation, n_total, decile_target, .drop=F) %>% tally(name="n_decile") %>%
  mutate(prop_decile=n_decile/n_total)
head(t2)
dim(t1)


t3 = bind_rows(mutate(t1, corr="mask"), mutate(t2, corr="no-mask")) %>%
  mutate(decile_target=factor(decile_target, levels=c(1:10), labels=paste0("dec",1:10)),
         corr=factor(corr, levels=c("mask","no-mask"), labels=c("mask dropouts", "no mask dropouts")))

ggplot(t3, aes(x=decile_target, y=prop_decile))+
  ggbeeswarm::geom_quasirandom(width = .2)+
  geom_boxplot(alpha=.5, color="red3", width=.5, outliers=F)+
  facet_wrap(vars(corr), ncol=1)+
  labs(title="Prop. module targets in expr. deciles",
       subtitle="13844 no low UMI logcounts")


t1 = filter(lg.mask, abs(rho)>.03, importance>.5) %>%
  group_by(TF, regulation) %>% add_tally(name="n_total") %>%
  filter(n_total>=20, target %in% mbv.degs) %>% 
  group_by(target, regulation) %>% tally(name="n_modules") %>%
  left_join(select(avg.expr, gene_name, "decile_target"=decile), by=c("target"="gene_name"))

t2 = filter(lg.nomask, abs(rho)>.03, importance>.5) %>%
  group_by(TF, regulation) %>% add_tally(name="n_total") %>%
  filter(n_total>=20, target %in% mbv.degs) %>% 
  group_by(target, regulation) %>% tally(name="n_modules") %>%
  left_join(select(avg.expr, gene_name, "decile_target"=decile), by=c("target"="gene_name"))

t3 = bind_rows(mutate(t1, corr="mask"), mutate(t2, corr="no-mask")) %>%
  mutate(decile_target=factor(decile_target, levels=c(1:10), labels=paste0("dec",1:10)),
         corr=factor(corr, levels=c("mask","no-mask"), labels=c("mask dropouts", "no mask dropouts")))

nrow(t1) #masked modules include every single F test gene (649 present)
nrow(t2) #modules without masking exclude some F test genes (573 present)

ggplot(t3, aes(x=decile_target, y=n_modules))+
  ggbeeswarm::geom_quasirandom(width = .2)+
  geom_boxplot(alpha=.5, color="red3", width=.5, outliers=F)+
  facet_wrap(vars(corr), ncol=1)+
  labs(title="# of modules for MBv F-test targets",
       subtitle="13844 no low UMI logcounts")

length(intersect(t1$target, mbv.degs2)) #masked modules include every single MBv DEG (387 present)
length(intersect(t2$target, mbv.degs2)) #modules without masking exclude ~50 MBv DEGs (323 present)

ggplot(filter(t3, target %in% mbv.degs2), aes(x=decile_target, y=n_modules))+
  ggbeeswarm::geom_quasirandom(width = .2)+
  geom_boxplot(alpha=.5, color="red3", width=.5, outliers=F)+
  facet_wrap(vars(corr), ncol=1)+
  labs(title="# of modules for MBv DEG targets",
       subtitle="13844 no low UMI logcounts")
