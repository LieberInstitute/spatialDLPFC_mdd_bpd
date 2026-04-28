setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(dplyr)
	library(ggplot2)
})

cpList = readRDS("plots/colorPalettes.rds")
load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")

colData(spe_pseudo)[["CRH"]] = logcounts(spe_pseudo)[rowData(spe_pseudo)$gene_name=="CRH",]

summ.df = as.data.frame(colData(spe_pseudo)[,c("condition", "sex", "smoothed_k9_1663", "CRH")])

df1 = bind_rows(filter(summ.df, #key_genes=="CRH", 
                       smoothed_k9_1663=="L1", condition=="NTC", sex=="F") %>%
                  mutate(cluster="c1", condition="ref1"),
                filter(summ.df, #key_genes=="CRH", 
                       smoothed_k9_1663=="L5", condition=="NTC", sex=="F") %>%
                  mutate(cluster="c2", condition="ref1"),
                filter(summ.df, #key_genes=="CRH", 
                       smoothed_k9_1663=="L6", condition=="NTC", sex=="F") %>%
                  mutate(cluster="c3", condition="ref1") %>%
                  filter(CRH>3)
)



df2 = bind_rows(filter(summ.df, #key_genes=="CRH", 
                       smoothed_k9_1663=="L2", condition=="NTC", sex=="F") %>%
                  mutate(cluster="c1", condition="test1", CRH= CRH-.4),
                filter(summ.df, #key_genes=="CRH",
                       smoothed_k9_1663=="L6", condition=="MDD", sex=="F") %>%
                  mutate(cluster="c2", condition="test1"),
                filter(summ.df, #key_genes=="CRH", 
                       smoothed_k9_1663=="L6", condition=="NTC", sex=="M") %>%
                  mutate(cluster="c3", condition="test1", CRH= CRH-.4) %>%
                  filter(CRH<5.65)
)

df3 = bind_rows(df1, df2) %>% 
  mutate(cluster= factor(cluster, levels=c("c1","c2","c3")),
         condition= factor(condition, levels=c("ref1","test1"), labels=c("reference","test")))

df4 = bind_rows(mutate(df3, cluster2= cluster), 
                mutate(df3, cluster2="all")) %>%
  mutate(cluster2= factor(cluster2, levels=c("all","c1","c2","c3")))

set.seed(123)
sub_rows = sample(1:nrow(df4), size=.8*nrow(df4))

df4 = mutate(df4[sub_rows,], CRH= CRH-.9)

# remove "all" and replace will re-aggregate
df4 = filter(df4, cluster2!="all")
df4 = bind_rows(df4, 
	mutate(df4, cluster2=factor("all", levels=c("all","c1","c2","c3")))
)

#just use this one
p1 <- ggplot(df4, aes(x=condition, y=CRH, shape=cluster))+
  ggbeeswarm::geom_quasirandom(size=2)+
  scale_shape_manual(values=c(1,0,2), guide="none")+
  scale_x_discrete(labels=c("NTC","DX"))+
  scale_y_continuous("logcounts", limits=c(1,8), breaks=c(1,3,5,7), labels=c("0","2","4","6"))+facet_grid(cols=vars(cluster2))+
  theme_bw()+theme(aspect.ratio=1, panel.grid.minor=element_blank())

ggsave("plots/publication/Figure1/LA-LR-model_fake-data.pdf",
       p1, width=5, height=3)


cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
