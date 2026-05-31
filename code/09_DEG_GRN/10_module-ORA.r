setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(clusterProfiler)
	library(ggplot2)
})

set.seed(123)

refined.modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv")

mod_subset = c("EEF1A1",#"EIF1",
        "PLP1","CD74","COX4I1","GLUL","IFITM3","GRIN1",
        "SNHG14","CAMK2N1","GAD1","A2M","PRKAR1A", #"UQCRH",
        #"ADAMTS1",
	"HSPA1A","MT1M",#"JUNB",
        "GFAP",
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

#stop("Early stopping to pick terms")

select.terms = list("A2M"=c("GO:0150104","GO:0003018","GO:2001214"), #A2M
                 "IFITM3"=c("GO:0071356","GO:0002526","GO:0019955"), #IFITM3
                 "CD74"=c("GO:0098883","GO:0001774","GO:0001818"), #CD74
                 "HSPA1A"=c("GO:2001233","GO:0048545","GO:0034599"), #HSPA1A "GO:1900745",
                 "MT1M"=c("GO:0010038","GO:0098754","GO:0045089"), #MT1M
                 "SNHG14"=c("GO:0001217","GO:0042393","GO:0141108"), #SNHG14
                 "GLUL"=c("GO:0051384","GO:0036293","GO:0070374"), #GLUL "GO:0071385"
                 "GAD1"=c("GO:0098982","GO:0033555","GO:0001664"), #GAD1
                 "CAMK2N1"=c("GO:0008331"), #CAMK2N1
                 "GRIN1"=c("GO:0048167","GO:0042752","GO:0014069"), #GRIN1
                 "PRKAR1A"=c("GO:0004860","GO:0031625"), #PRKAR1A
                 "COX4I1"=c("GO:0046034","GO:0004129","GO:0032543"),#COX4I1
                 #"GO:0006119","GO:0007005","GO:0006979","GO:0006091",#COX4I1
                 "EEF1A1"=c("GO:0002181","GO:0043484"),#EEF1A1
                 "GFAP"=c("GO:0072331","GO:0010506","GO:1904018"), #GFAP
                 "PLP1"=c("GO:0008366","GO:0043209"), #PLP1
                 # DEGs missing from modules
                 "missing"=c("GO:0005179","GO:0004725","GO:0005833","GO:0062023","GO:0101002")
)

mod_subset = c("SNHG14","PRKAR1A","COX4I1","EEF1A1",
  "PLP1",
  "GFAP","GLUL","HSPA1A","MT1M",
  "IFITM3","CD74","A2M",
  "GRIN1","CAMK2N1","GAD1",
  "missing")
names(mod_subset) = mod_subset

#everyone rep 3 except GRIN1 and UQCRH
facet_y = do.call(c, lapply(mod_subset,
                            function(x) {
                              if(x %in% c("PRKAR1A","EEF1A1","PLP1")) {
                                return(rep(x,2))
                              } else {
                                if(x %in% c("CAMK2N1")) {
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
  tmp1 = tmp1[intersect(unlist(select.terms), rownames(tmp1)),] %>% mutate(module=x)
  return(tmp1)
}))


plot.df = mutate(plot.df, ID= factor(ID, levels=rev(unlist(select.terms[mod_subset]))), 
                 module= factor(module, levels=mod_subset),
                 #                 module2= factor(module, levels=rev(names(mod_subset))),
                 facet_rows= factor(ID, levels=unlist(select.terms[mod_subset]), labels=facet_y)) %>%
  arrange(ID)

plot.df$Description2 = factor(plot.df$Description, levels=unique(plot.df$Description))

p2 <- ggplot(plot.df, aes(y=Description2, x=module, size=Count, color=log2(FoldEnrichment)))+
  geom_count()+
  scale_color_gradient("log2\nFold\nEnrich.", low="white", high="black", limits=c(0,8))+
  scale_size("# DEGs", range=c(2,6), breaks=c(3,9,15,21))+
  facet_grid(rows=vars(facet_rows), scales="free_y", space="free_y")+
  labs(y="")+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
                   strip.text.y=element_text(angle=0))

ggsave(file="plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_subset-refined-modules_ORA-GO-dotplot.pdf",
	p2, width=7, height=10)
cat("\nSaved representative term dotplot to: plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_subset-refined-modules_ORA-GO-dotplot.pdf\n")


# how consistent is the GO phenotype with the DEGs in the module
resList = lapply(names(ora_modules), function(x) {
  res1 = ora_modules[[x]]@result
  t1 = strsplit(res1$geneID, "/")
  if(length(t1)>1) {
    t2 = unique(do.call(c, t1))
  } else {
    t2 = unique(t1[[1]])
  }
  return(t2)
})
names(resList) <- names(ora_modules)

totalList = lapply(names(ora_modules), function(x) {
  if(x=="missing") {
    return(length(missing.degs))
  } else {
    return(nrow(filter(refined.modules, TF==x))+1)
  }
})
names(totalList) = names(ora_modules)

df1 = data.frame("module"=names(ora_modules), "n_genes"=unlist(totalList), 
           "prop_GO"=sapply(names(ora_modules), function(x) round(length(resList[[x]])/totalList[[x]],2)))

df1$mod_group = factor(df1$module, levels=names(mod_subset),
                       labels=c("multi","multi","multi","multi",
				"Oligo",
                                "Astro","Astro",
				"Vasc","Vasc","Vasc","Micro","Vasc",
                                "ExcN","ExcN","InhN",
                                "multi"))

col.pal = c("Vasc"=cpList$low.res.light[["Micro.Vasc"]],
            "Micro"=cpList$low.res.light[["L3"]],
            cpList$low.res.light[c("Astro","Oligo")],
            "InhN"=cpList$low.res.light[["Inhb"]],
            "ExcN"=cpList$low.res.light[["L2"]],
            "multi"="grey"
)

p1 <- ggplot(df1, aes(x=n_genes, y=prop_GO, fill=mod_group))+
  ggrepel::geom_label_repel(aes(label=module), min.segment.length=0)+
  geom_point(shape=21, size=2)+scale_fill_manual(values=col.pal)+
  #geom_label(aes(label=module), nudge_y = .015)+
  coord_cartesian(xlim=c(0,120), ylim=c(0,1))+
  labs(x="Module size (# DEGs)", y="Prop. of module represented\nin sig. GO terms")+
  theme_bw()+theme(panel.grid.minor=element_blank(), aspect.ratio=1,
                   legend.position="none")

ggsave(file="plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_subset-refined-modules_ORA-GO-scatter.pdf",
       p1, width=5, height=5)

## Reproducibility information
cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
