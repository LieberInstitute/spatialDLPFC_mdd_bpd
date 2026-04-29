setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(clusterProfiler)
	library(ggplot2)
})

set.seed(123)

refined.modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv")

mod_subset = c("A2M","IFITM3","CD74","HSPA1A","MT1X","SNHG14","GLUL",
               "CAMK2N1","GAD1","GRIN1","PRKAR1A","UQCRH",
               "APLP1","FTL","PLP1",
		"missing")
names(mod_subset) = mod_subset

# DEGs missing from modules
source("code/09_DEG_GRN/load_DEGs.r")
all.degs = unique(do.call(rbind, sigList)$gene_name)
cat("\nNumber of DEGs total:", length(all.degs),"\n")
#should be 503
missing.degs = setdiff(all.degs, union(refined.modules$TF, refined.modules$target))
cat("Number of DEGs missing from modules:", length(missing.degs),"\n")
#should be 66

if(file.exists("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_modules-DEG-subset-refined_ORA-GO-results.rda")) {
	ora_modules = readRDS("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_modules-DEG-subset-refined_ORA-GO-results.rda")
	cat("\nLoaded saved ORA GO results from: processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_modules-DEG-subset-refined_ORA-GO-results.rda\n")
} else {

	ora_modules = lapply(mod_subset, function(x) {
	  set.seed(123)
	  if(x=="missing") {
		mod.genes=missing.degs
	  } else {mod.genes = c(filter(refined.modules, TF==x)$target, x)}
	  ora = enrichGO(mod.genes, "org.Hs.eg.db", keyType="SYMBOL", ont="ALL", pAdjustMethod="BH")#,
	                 #universe=gene_universe)
	  return(ora)
	})

	cat("\nNumber of sig. terms per module:\n")
	sapply(ora_modules, dim)

	saveRDS(ora_modules, "processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_modules-DEG-subset-refined_ORA-GO-results.rda")
	cat("\nSaved ORA GO results to: processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_modules-DEG-subset-refined_ORA-GO-results.rda\n")
}

# plot select terms
select.terms = c("GO:0150104","GO:0003018","GO:2001214", #A2M
                 "GO:0071356","GO:0002526","GO:0035456", #IFITM3
                 "GO:0098883","GO:0042116","GO:0002443", #CD74
                 "GO:0034605","GO:0048545","GO:0034599", #HSPA1A "GO:1900745",
                 "GO:0010038","GO:0050727","GO:0098754", #MT1X
		 "GO:0051010", #SNHG14
                 "GO:0070371","GO:0003158","GO:0001666", #GLUL "GO:0071385"
		 "GO:0008331", #CAMK2N1
                 "GO:0098982","GO:0005184", #GAD1 #"GO:0060077"
                 "GO:0014069", "GO:0034399", #GRIN1
                 "GO:0031625","GO:0019207",#PRKAR1A #"GO:0034236", "GO:0004860",
                 "GO:0045333","GO:0007005", #UQCRH
                 "GO:0031109","GO:0008380","GO:0031690", #APLP1
                 "GO:0002181","GO:0015934","GO:0015935", #FTL
                 #"GO:0005840","GO:0030863","GO:0005925", #FTL
                 "GO:0008366","GO:0007163","GO:0032535", #PLP1
		 # DEGs missing from modules
		 "GO:0005179","GO:0004725","GO:0005833","GO:0062023","GO:0101002"
                 )

#everyone rep 3 except GRIN1 and UQCRH
facet_y = do.call(c, lapply(mod_subset,
       function(x) {
         if(x %in% c("GAD1","GRIN1","PRKAR1A","UQCRH")) {
           return(rep(x,2))
         } else {
	   if(x %in% c("SNHG14","CAMK2N1")) {
	      return(x)
	   } else {
		if(x=="missing") {
		   return(rep(x,5))
		} else {
		   return(rep(x,3))
		}
	   }
         }
       })
)

plot.df = do.call(rbind, lapply(mod_subset, function(x) {
  tmp1 = ora_modules[[x]]@result
  tmp1 = tmp1[intersect(select.terms, rownames(tmp1)),] %>% mutate(module=x)
  return(tmp1)
}))

plot.df = mutate(plot.df, ID= factor(ID, levels=rev(select.terms)), 
                 module= factor(module, levels=names(mod_subset)),
#                 module2= factor(module, levels=rev(names(mod_subset))),
                 facet_rows= factor(ID, levels=select.terms, labels=facet_y)) %>%
  arrange(ID)

plot.df$Description2 = factor(plot.df$Description, levels=unique(plot.df$Description))

p1 <- ggplot(plot.df, aes(y=Description2, x=module, size=Count, color=log2(FoldEnrichment)))+
  geom_count()+
  scale_color_gradient("log2\nFold\nEnrich.", low="white", high="black", limits=c(0,8))+
  scale_size("# DEGs", range=c(2,6), breaks=c(3,9,15,21))+
  facet_grid(rows=vars(facet_rows), scales="free_y", space="free_y")+
  labs(y="")+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
		strip.text.y=element_text(angle=0))

ggsave(file="plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_subset-refined-modules_ORA-GO-dotplot.pdf",
	p1, width=7, height=10)
cat("\nSaved representative term dotplot to: plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_subset-refined-modules_ORA-GO-dotplot.pdf\n")


## Reproducibility information
cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
