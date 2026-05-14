setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
})

set.seed(123)

regulons <- read.csv("processed-data/10_SCENIC/spe-n109_13162-no-lowUMI_logcounts_regulons-top20_merged.csv", row.names=1) 

#check1 = filter(regulons, set_size>5)$TF

reg_subset = c("KLF4","JUNB","FOS","FOSL2","CEBPD",
               "JUN","SOX2","SOX9","ETS2","THRA",
               "NR2F2","DLX1","ARX","LHX6","TFEB",
               "SOX8","SOX10","NFIA","SREBF1","BCL6")
#setdiff(check1, reg_subset) #they are the same
names(reg_subset) = reg_subset
reg_subset = as.list(reg_subset)


source("code/09_DEG_GRN/load_DEGs.r")

all.degs = unique(sig.df$gene_name)

sigList[["sm.la"]]$cluster = as.factor(sigList[["sm.la"]]$cluster)
sigList[["se.la"]]$cluster = as.factor(sigList[["se.la"]]$cluster)

regList <- lapply(reg_subset, function(y) {
  reg.genes = c(unlist(strsplit(filter(regulons, TF==y)$set_str, "/")), y)
  t1 = do.call(rbind, lapply(sigList[c("sm.la","sm.lr")], function(x) {
    filter(x, gene_name %in% reg.genes) %>% 
      distinct(dir, cluster, gene_name, gene_id) %>%
      group_by(dir, cluster, .drop=F) %>% tally()
    })) %>%
    mutate(annot="sm")
  t2 = do.call(rbind, lapply(sigList[c("se.la","se.lr")], function(x) {
    filter(x, gene_name %in% reg.genes) %>%
      distinct(dir, cluster, gene_name, gene_id) %>%
      group_by(dir, cluster, .drop=F) %>% tally()
    })) %>%
    mutate(annot="se")
  
  bind_rows(t1, t2) %>% mutate(regulon=y, prop_reg=n/(length(reg.genes)-1)) 
})

clus.df = do.call(rbind, regList) %>%
  mutate(is_LA= factor(ifelse(cluster=="L-A", "L-A","L-R"), levels=rev(c("L-A","L-R"))), 
         annot= factor(annot, levels=c("sm","se")),
         cluster= factor(cluster, levels=c("L-A","Micro.Vasc","Astro","L1","L2","L2.3","L3.4","L4",
                                          "Inhb","L5","L6","WM","Oligo"),
                         labels=c("L-A","M.V","Ast","L1","L2","L2.3","L3.4","L4",
                                  "Inb","L5","L6","WM","Olg")),
         regulon= factor(regulon, levels=names(reg_subset)),
         sign_n= ifelse(dir=="Dec.", -n, n),
         sign_prop= ifelse(dir=="Dec.", -prop_reg, prop_reg))

col.pal = rev(RColorBrewer::brewer.pal("RdBu", n=8))
col.pal = c(col.pal[1:4], "white", col.pal[5:8])

p1 <- ggplot(clus.df, aes(x=cluster, y=dir, fill=sign_prop))+
  geom_tile(color="black", linewidth=.1)+
  geom_text(data=filter(clus.df, n>2), aes(label=n), size=2)+
  scale_fill_gradientn(colors=col.pal, limits=c(-1,1))+
  facet_grid(rows=vars(regulon), cols=vars(annot), scales="free_x", space="free_x")+
  labs(fill="prop.\nregulon\nDEGs", x="", y="regulon")+
  theme_minimal()+theme(strip.text.y=element_text(angle=0, hjust=0))

# by dx*sex group
regList2 <- lapply(reg_subset, function(y) {
  reg.genes = c(unlist(strsplit(filter(regulons, TF==y)$set_str, "/")), y)
  t1 = filter(sigList[["sm.la"]], gene_name %in% reg.genes) %>%
    group_by(dir, sex.group, .drop=F) %>% tally() %>%
    mutate(is_LA="L-A", annot="sm")
  t2 = filter(sigList[["sm.lr"]], gene_name %in% reg.genes) %>%
    group_by(dir, sex.group, gene_name, .drop=F) %>% tally() %>%
    filter(!is.na(gene_name)) %>%
    group_by(dir, sex.group, .drop=F) %>% tally() %>%
    mutate(is_LA="L-R", annot="sm")
  
  t3 = filter(sigList[["se.la"]], gene_name %in% reg.genes) %>%
    group_by(dir, sex.group, .drop=F) %>% tally() %>%
    mutate(is_LA="L-A", annot="se")
  t4 = filter(sigList[["se.lr"]], gene_name %in% reg.genes) %>%
    group_by(dir, sex.group, gene_name, .drop=F) %>% tally() %>%
    filter(!is.na(gene_name)) %>%
    group_by(dir, sex.group, .drop=F) %>% tally() %>%
    mutate(is_LA="L-R", annot="se")
  
  bind_rows(t1, t2, t3, t4) %>% mutate(regulon=y, prop_reg=n/(length(reg.genes)-1)) 
})

dxsex.df = do.call(rbind, regList2) %>%
  mutate(annot= factor(annot, levels=c("sm","se")),
         is_LA= factor(is_LA, levels=rev(c("L-A","L-R"))),
         regulon= factor(regulon, levels=names(reg_subset)),
         facet_col= factor(paste(dir, annot), levels=c("Dec. sm","Dec. se","Inc. sm","Inc. se")),
         sign_n= ifelse(dir=="Dec.", -n, n),
         sign_prop= ifelse(dir=="Dec.", -prop_reg, prop_reg))


p2 <- ggplot(dxsex.df, aes(x=facet_col, y=is_LA, fill=sign_prop))+
  geom_tile(color="black", linewidth=.1)+
  geom_text(data=filter(dxsex.df, n>2), aes(label=n), size=2)+
  scale_fill_gradientn(colors=col.pal, limits=c(-1,1))+
  facet_grid(rows=vars(regulon), cols=vars(sex.group))+#, scales="free_x", space="free_x")+
  labs(fill="prop.\nregulon\nDEGs", x="", y="DE model")+
  theme_minimal()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
	strip.text.y=element_text(angle=0, hjust=0))

pdf(file="plots/10_SCENIC/spe-n109_13162-no-lowUMI_regulons-top20_DE-group-plot.pdf",
	height=6, width=6)
p1
p2
dev.off()
cat("\nPlots saved to: plots/10_SCENIC/spe-n109_13162-no-lowUMI_regulons-top20_DE-group-plot.pdf\n")



## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
