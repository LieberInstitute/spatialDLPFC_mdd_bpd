setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(ggrastr)
})

source("code/publication/plotting_utils.r")
length(unique(sig.df$gene_id)) #503

mratio.sn = read.csv("processed-data/06_pseudobulk/SZBDMulti-seq/SZBD-control_azimuth-super-broad_mean-ratio.csv") %>%
  group_by(gene_name) %>% slice_max(MeanRatio, n=1) %>%
  mutate(cellType.target= ifelse(MeanRatio<1.5, "multi", cellType.target)) %>% ungroup()

df = mutate(as.data.frame(colData(sce_summ)), 
            azimuth_super.broad=factor(azimuth_super.broad, 
                                       levels=c("Vasc","Astro","Oligo","Micro","InhN","ExcN")))

plot.genes = c("ITIH5","COL5A3","PLP1","C3","GAD1","RALYL","SLC6A1")

df1 = do.call(rbind, lapply(plot.genes, function(x) {
  tmp = df
  tmp$gene_name = x
  tmp$avg.logcounts = assay(sce_summ, "logcounts.mean")[rowData(sce_summ)$gene_name==x,]
  tmp$prop.detected = assay(sce_summ, "logcounts.prop.detected")[rowData(sce_summ)$gene_name==x,]
  return(tmp)
})) %>% mutate(gene_name= factor(gene_name, levels=plot.genes))

plist1 <- lapply(plot.genes, function(i) {
  ymax1 = as.numeric(as.character(factor(i, levels=plot.genes, labels=c(4,4,6,4,4,6,4))))
  ggplot(filter(df1, gene_name==i), aes(x=azimuth_super.broad, y=avg.logcounts, fill=azimuth_super.broad))+
    geom_bar(stat="identity", color="black", linewidth=.1)+scale_fill_manual(values=sn.col.pal, guide="none")+
    labs(title=i, y="mean expr.")+scale_y_continuous(limits=c(0,ymax1), breaks=seq(0, ymax1, length.out=3))+
    theme_minimal()+theme(text=element_text(size=6), plot.title=element_text(size=7), axis.text.x=element_text(angle=90, hjust=1, vjust=.5), 
                     axis.title.x=element_blank(), panel.grid.major.x=element_blank())
})

ggsave(file="plots/publication/supp_mean-ratio/expr-logcounts_bar.pdf", 
       marrangeGrob(plist1, ncol=1, nrow=1, top=NULL),
       height=1, width=1)




mbar <- getMeanRatioBar(plot.genes, sce_summ)+theme(axis.text.y=element_text(size=8))
ggsave(file="plots/publication/supp_mean-ratio/expr-fill_bar.pdf", mbar, height=3.5, width=3)




plot.df = mutate(sig.df, cluster_source= factor(cluster_source, levels=all_clusters),
                 consensus_group= factor(cluster_source, levels=all_clusters, 
                                         labels=c("L-A","L-A","M.V","Ast","L1",
                                                  #"L2","L2.3","L3.4","L4",
                                                  "L2/3","L2/3","L3/4","L3/4",
                                                  "Inb",
                                                  "L5","L5","L6","L6","WM/O","WM/O")),
                 source=factor(source, levels=c("sm","se"), labels=c("domain-SP","domain-CT")))

#total.degs = distinct(plot.df, source, cluster_source, gene_name) %>% group_by(source, cluster_source) %>% tally()

plot.df = left_join(plot.df, mratio.sn[,c("gene_name","gene","MeanRatio","cellType.target")], by=c("gene_name","gene_id"="gene")) %>%
  filter(!is.na(cellType.target)) %>%
  mutate(cellType.target=factor(cellType.target, levels=c("multi","ExcN","InhN","Micro","Oligo","Astro","Vasc")))

#sn.degs = distinct(plot.df, source, cluster_source, gene_name) %>% group_by(source, cluster_source) %>% tally()

tally.df = group_by(plot.df, consensus_group, source, cluster_source, cellType.target) %>% tally()

#use darker oligo color palette because of white outline (consistent with Fig2)
color.palette = c("Vasc"=cpList$low.res.light[["Micro.Vasc"]],
            "Micro"="#C28658",
            cpList$low.res.light[c("Astro")], cpList$low.res.bright["Oligo"],
            "InhN"=cpList$low.res.light[["Inhb"]],
            "ExcN"=cpList$low.res.light[["L2"]],
            "multi"="grey"
)

p1 <- ggplot(tally.df, aes(x=cluster_source, y=n, fill=cellType.target))+
  geom_bar(stat="identity", position="fill", color="white", linewidth=.3)+
  scale_fill_manual(values=color.palette)+
  facet_wrap(vars(consensus_group), scales="free_x")+
  theme_minimal()+theme(text=element_text(size=6), axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
                        legend.position="none")

ggsave(file="plots/publication/supp_mean-ratio/mratio-tally_consensus-wrap.pdf", p1, height=4, width=3)

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
