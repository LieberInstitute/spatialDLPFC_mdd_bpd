setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
})

set.seed(123)

refined.modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv")

#mod_subset = c("A2M","IFITM3","CD74","HSPA1A","MT1X","SNHG14","GLUL",
#               "CAMK2N1","GAD1","GRIN1","PRKAR1A","UQCRH",
#               "APLP1","FTL","PLP1")
mod_subset = c("SNHG14","PRKAR1A","COX4I1","EEF1A1",
  "PLP1",
  "GFAP","GLUL","HSPA1A","MT1M",
  "IFITM3","CD74","A2M",
  "GRIN1","CAMK2N1","GAD1")
names(mod_subset) = mod_subset
mod_subset = as.list(mod_subset)


source("code/09_DEG_GRN/load_DEGs.r")

all.degs = unique(do.call(rbind, sigList)$gene_name)
missing.degs = setdiff(all.degs, union(refined.modules$TF, refined.modules$target))

sigList[["sm.la"]]$cluster = as.factor(sigList[["sm.la"]]$cluster)
sigList[["se.la"]]$cluster = as.factor(sigList[["se.la"]]$cluster)

mod_subset[["missing"]] = missing.degs

modList <- lapply(mod_subset, function(y) {
  if(length(y)==1) {
	mod.genes = c(filter(refined.modules, TF==y)$target, y)
  } else {
	mod.genes = y
	y = "missing"
  }
  t1 = do.call(rbind, lapply(sigList[c("sm.la","sm.lr")], function(x) {
    filter(x, gene_name %in% mod.genes) %>% 
      distinct(dir, cluster, gene_name, gene_id) %>%
      group_by(dir, cluster, .drop=F) %>% tally()
    })) %>%
    mutate(annot="sm")
  t2 = do.call(rbind, lapply(sigList[c("se.la","se.lr")], function(x) {
    filter(x, gene_name %in% mod.genes) %>%
      distinct(dir, cluster, gene_name, gene_id) %>%
      group_by(dir, cluster, .drop=F) %>% tally()
    })) %>%
    mutate(annot="se")
  
  bind_rows(t1, t2) %>% mutate(module=y, prop_mod=n/length(mod.genes)) 
})

clus.df = do.call(rbind, modList) %>%
  mutate(is_LA= factor(ifelse(cluster=="L-A", "L-A","L-R"), levels=rev(c("L-A","L-R"))), 
         annot= factor(annot, levels=c("sm","se")),
         cluster= factor(cluster, levels=c("L-A","Micro.Vasc","Astro","L1","L2","L2.3","L3.4","L4",
                                          "Inhb","L5","L6","WM","Oligo"),
                         labels=c("L-A","M.V","Ast","L1","L2","L2.3","L3.4","L4",
                                  "Inb","L5","L6","WM","Olg")),
         module= factor(module, levels=names(mod_subset)),
         sign_n= ifelse(dir=="Dec.", -n, n),
         sign_prop= ifelse(dir=="Dec.", -prop_mod, prop_mod))

col.pal = rev(RColorBrewer::brewer.pal("RdBu", n=8))
col.pal = c(col.pal[1:4], "white", col.pal[5:8])

p1 <- ggplot(clus.df, aes(x=cluster, y=dir, fill=sign_prop))+
  geom_tile(color="black", linewidth=.1)+
  geom_text(data=filter(clus.df, n>2), aes(label=n), size=2)+
  scale_fill_gradientn(colors=col.pal, limits=c(-1,1))+
  facet_grid(rows=vars(module), cols=vars(annot), scales="free_x", space="free_x")+
  labs(fill="prop.\nmodule\nDEGs", x="", y="DEG module")+
  theme_minimal()+theme(strip.text.y=element_text(angle=0, hjust=0))

# by dx*sex group
modList2 <- lapply(mod_subset, function(y) {
    if(length(y)==1) {
        mod.genes = c(filter(refined.modules, TF==y)$target, y)
  } else {
        mod.genes = y
        y = "missing"
  }
  t1 = filter(sigList[["sm.la"]], gene_name %in% mod.genes) %>%
    group_by(dir, sex.group, .drop=F) %>% tally() %>%
    mutate(is_LA="L-A", annot="sm")
  t2 = filter(sigList[["sm.lr"]], gene_name %in% mod.genes) %>%
    group_by(dir, sex.group, gene_name, .drop=F) %>% tally() %>%
    filter(!is.na(gene_name)) %>%
    group_by(dir, sex.group, .drop=F) %>% tally() %>%
    mutate(is_LA="L-R", annot="sm")
  
  t3 = filter(sigList[["se.la"]], gene_name %in% mod.genes) %>%
    group_by(dir, sex.group, .drop=F) %>% tally() %>%
    mutate(is_LA="L-A", annot="se")
  t4 = filter(sigList[["se.lr"]], gene_name %in% mod.genes) %>%
    group_by(dir, sex.group, gene_name, .drop=F) %>% tally() %>%
    filter(!is.na(gene_name)) %>%
    group_by(dir, sex.group, .drop=F) %>% tally() %>%
    mutate(is_LA="L-R", annot="se")
  
  bind_rows(t1, t2, t3, t4) %>% mutate(module=y, prop_mod=n/length(mod.genes)) 
})

dxsex.df = do.call(rbind, modList2) %>%
  mutate(annot= factor(annot, levels=c("sm","se")),
         is_LA= factor(is_LA, levels=rev(c("L-A","L-R"))),
         module= factor(module, levels=names(mod_subset)),
         facet_col= factor(paste(dir, annot), levels=c("Dec. sm","Dec. se","Inc. sm","Inc. se")),
         sign_n= ifelse(dir=="Dec.", -n, n),
         sign_prop= ifelse(dir=="Dec.", -prop_mod, prop_mod))


p2 <- ggplot(dxsex.df, aes(x=facet_col, y=is_LA, fill=sign_prop))+
  geom_tile(color="black", linewidth=.1)+
  geom_text(data=filter(dxsex.df, n>2), aes(label=n), size=2)+
  scale_fill_gradientn(colors=col.pal, limits=c(-1,1))+
  facet_grid(rows=vars(module), cols=vars(sex.group))+#, scales="free_x", space="free_x")+
  labs(fill="prop.\nmodule\nDEGs", x="", y="DE model")+
  theme_minimal()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
	strip.text.y=element_text(angle=0, hjust=0))

pdf(file="plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_subset-refined-modules_DE-group-plot.pdf",
	height=6, width=6)
p1
p2
dev.off()
cat("\nPlots saved to: plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_subset-refined-modules_DE-group-plot.pdf\n")



## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
