setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(clusterProfiler)
	library(ggplot2)
})

set.seed(123)

ora_modules = readRDS("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_modules-DEG-subset-refined_ORA-GO-results.rda")

select.terms = list("A2M"=c("GO:0150104","GO:0003018","GO:2001214"), #A2M
                    "IFITM3"=c("GO:0062023","GO:0002526","GO:0019955"), #IFITM3
                    "CD74"=c("GO:0098883","GO:0001774","GO:0001818"), #CD74
                    "HSPA1A"=c("GO:0048545","GO:0034599","GO:0034612"), #HSPA1A
                    "MT1M"=c("GO:0010038","GO:0098754","GO:0045089"), #MT1M
                    "SNHG14"=c("GO:0001217","GO:0042393","GO:0141108"), #SNHG14
                    "GLUL"=c("GO:0036293","GO:0070374","GO:0051384"), #GLUL
                    "GAD1"=c("GO:0098982","GO:0033555","GO:0160041"), #GAD1
                    "CAMK2N1"=c("GO:0008331"), #CAMK2N1
                    "GRIN1"=c("GO:0048167","GO:0042752","GO:0014069"), #GRIN1
                    "PRKAR1A"=c("GO:0004860","GO:0031625"), #PRKAR1A
                    "COX4I1"=c("GO:0046034","GO:0004129","GO:0032543"),#COX4I1
                    "EEF1A1"=c("GO:0002181","GO:0043484"),#EEF1A1
                    "GFAP"=c("GO:0072331","GO:0010506","GO:1904018"), #GFAP
                    "PLP1"=c("GO:0008366","GO:0043209")#, #PLP1
)


mod_subset = c("SNHG14","PRKAR1A","COX4I1","EEF1A1",
  "PLP1",
  "GFAP","GLUL","HSPA1A","MT1M",
  "IFITM3","CD74","A2M",
  "GRIN1","CAMK2N1","GAD1")
names(mod_subset) = mod_subset

facet_y = do.call(c, lapply(mod_subset,
                            function(x) {
                              if(x %in% c("PRKAR1A","EEF1A1","PLP1")) {
                                return(rep(x,2))
                              } else {
                                if(x %in% c("CAMK2N1")) {
                                  return(x)
                                } else {
                                    return(rep(x,3))
                                  }
                              }
                            })
)

plot.df = do.call(rbind, lapply(mod_subset, function(x) {
  tmp1 = ora_modules[[x]]@result
  tmp1 = tmp1[intersect(unlist(select.terms), rownames(tmp1)),] %>% mutate(module=x)
  return(tmp1)
}))


plot.df = mutate(plot.df, ID= factor(ID, levels=rev(unlist(select.terms[mod_subset]))), 
                 module= factor(module, levels=mod_subset),
                 facet_rows= factor(ID, levels=unlist(select.terms[mod_subset]), labels=facet_y)) %>%
  arrange(ID)

plot.df$Description2 = factor(plot.df$Description, levels=unique(plot.df$Description))

# print limits to output so can set the same for missing DEG ORA plot
summary(plot.df$Count)

p1 <- ggplot(plot.df, aes(y=Description2, x=module, size=Count, color=log2(FoldEnrichment)))+
  geom_count()+
  scale_color_gradient("log2\nFold\nEnrich.", low="white", high="black", limits=c(0,8))+
  scale_size("# DEGs", range=c(2,6), breaks=c(3,9,15), limits=c(2,19))+
  facet_grid(rows=vars(facet_rows), scales="free_y", space="free_y")+
  labs(y="")+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
                   strip.text.y=element_text(angle=0))

ggsave(file="plots/publication/Figure2/modules_GO-ORA_dotplot.pdf", p1, width=7, height=9)

## Reproducibility information
cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()

