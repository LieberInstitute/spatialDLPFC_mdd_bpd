setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
})

t1 = read.table("processed-data/00_genotypes/GWAS-stats/gwas_panel_metrics.tab", header=T)
t1$disorder = factor(t1$disorder, levels=c("MDD","BD","SCZD"), labels=c("MDD","BD","SCZ"))

p1 <- ggplot(t1.1, aes(x=disorder, y=count, fill=grouping))+
  geom_bar(stat="identity", position="dodge", width=.6)+
  geom_text(aes(label=count), color="red", hjust=.5, size=2)+
  scale_fill_manual(values=c("grey50","black"))+
  scale_y_continuous("European participants", labels=function(x) format(x, scientific=T))+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(), text=element_text(size=6),
                        legend.key.size=unit(6,"pt"))

ggsave(file="plots/publication/Figure_eQTL/supp_GWAS-participants.pdf", p1,
       width=3, height=1.5)

# venn diagrams created by hand using this file: processed-data/00_genotypes/GWAS-stats/gwas_variant_overlaps.tab

t2 = read.table("processed-data/00_genotypes/GWAS-stats/gwas_variant_overlaps.tab", header=T)

t2.1 = tidyr::pivot_longer(t2[,1:7], c("n_bd","n_mdd","n_sczd"), names_to="disorder", values_to="variants",
                           names_prefix="n_") %>%
  filter(scheme!="explore_all") %>%
  mutate(disorder= factor(disorder, levels=c("mdd","bd","sczd"), labels=c("MDD","BD","SCZ")),
         scheme= factor(scheme, levels=c("strict_all","mood_explore")))

p2 <- ggplot(t2.1, aes(x=disorder, y=variants))+
  geom_bar(stat="identity")+
  facet_wrap(vars(scheme))+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
                        text=element_text(size=6))

ggsave(file="plots/publication/Figure_eQTL/supp_GWAS-variants.pdf", p2,
       width=2, height=1)


## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
