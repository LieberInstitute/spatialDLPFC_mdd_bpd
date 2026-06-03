setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
})

set.seed(123)

refined.modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv")

mod_subset = c("SNHG14","PRKAR1A","COX4I1","EEF1A1",
  "PLP1",
  "GFAP","GLUL","HSPA1A","MT1M",
  "IFITM3","CD74","A2M",
  "GRIN1","CAMK2N1","GAD1")
names(mod_subset) = mod_subset
mod_subset = as.list(mod_subset)


source("code/09_DEG_GRN/load_DEGs.r")
mbv.degs = unique(sig.df$gene_name)
missing.degs = setdiff(mbv.degs, union(refined.modules$TF, refined.modules$target))

sigList[["sm.la"]]$cluster = as.factor(sigList[["sm.la"]]$cluster)
sigList[["se.la"]]$cluster = as.factor(sigList[["se.la"]]$cluster)

mod_subset[["missing"]] = missing.degs

# by domain
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

# by dx*sex 
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


clus.df = do.call(rbind, modList) %>%
  mutate(y_values= factor(dir, levels=c("Dec.","Inc.")),
         facet_col= factor(annot, levels=c("sm","se"), labels=c("domain-SP","domain-CT")),
         x_values= factor(cluster, levels=c("L-A","Micro.Vasc","Astro","L1","L2","L2.3","L3.4","L4",
                                           "Inhb","L5","L6","WM","Oligo"),
                         labels=c("L-A","M.V","Ast","L1","L2","L2.3","L3.4","L4",
                                  "Inb","L5","L6","WM","Olg")),
         module= factor(module, levels=names(mod_subset)),
         sign_n= ifelse(dir=="Dec.", -n, n),
         sign_prop= ifelse(dir=="Dec.", -prop_mod, prop_mod))

dxsex.df = do.call(rbind, modList2) %>%
  mutate(x_values= factor(paste(is_LA, annot), levels=c("L-A sm","L-A se","L-R sm","L-R se"),
                          labels=c("L-A domain-SP","L-A domain-CT","L-R domain-SP","L-R domain-CT")),
         y_values= factor(dir, levels=c("Dec.","Inc.")),
         module= factor(module, levels=names(mod_subset)),
         sign_n= ifelse(dir=="Dec.", -n, n),
         sign_prop= ifelse(dir=="Dec.", -prop_mod, prop_mod))

joined.df = bind_rows(select(clus.df, y_values, x_values, facet_row=module, facet_col, sign_prop, n),
          select(dxsex.df, y_values, x_values, facet_row=module, facet_col=sex.group, sign_prop, n)) %>%
  mutate(y_values= factor(y_values, levels=c("Dec.","Inc.")), 
         x_values= factor(x_values, levels=c(levels(clus.df$x_values), 
                                             "L-A domain-SP","L-A domain-CT","L-R domain-SP","L-R domain-CT")),
         facet_col= factor(facet_col, levels=c(levels(clus.df$facet_col), levels(dxsex.df$sex.group)),
                           labels=c(levels(clus.df$facet_col), gsub("_","\n",levels(dxsex.df$sex.group))))
  )

col.pal = rev(RColorBrewer::brewer.pal("RdBu", n=8))
col.pal = c(col.pal[1:4], "white", col.pal[5:8])

p1 <- ggplot(joined.df, aes(x=x_values, y=y_values, fill=sign_prop))+
  geom_tile(color="black", linewidth=.1)+
  geom_text(data=filter(joined.df, n>0), aes(label=n), size=2)+
  scale_fill_gradientn(colors=col.pal, limits=c(-1,1))+
  facet_grid(rows=vars(facet_row), cols=vars(facet_col), scales="free", space="free")+
  labs(fill="prop.\nmodule\nDEGs", y="sign(logFC)")+
  theme_minimal()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
                        strip.text.y=element_text(angle=0, hjust=0, size=8), 
                        strip.text.x=element_text(size=6), axis.title.x=element_blank())

pdf(file="plots/publication/modules/supp_module-DEG-groups.pdf",
    width=6.5, height=8)
p1+theme(legend.position="none")
p1
dev.off()

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()

