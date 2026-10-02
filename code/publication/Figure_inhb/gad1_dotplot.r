setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
})
set.seed(123)

refined.modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv")

source("code/publication/plotting_utils.r")
de.df2 = filter(de.df, group!="MDD.BPD")

mod_gene="GAD1"
#plot.genes = c(filter(refined.modules, TF==mod_gene)$target, mod_gene)

#ordered based on general vs lineage
plot.genes = c("PVALB","RGS5","KCNS3","TAC1","TRBC2",
"GAD1","GAD2","ZNF385D","PNOC","LGI2","SLC32A1","SLC6A1",
"DLX6-AS1","RELN","VIP")

#o1 = filter(refined.modules, TF==mod_gene, target %in% plot.genes) %>% arrange(desc(importance))
#o2 = filter(refined.modules, TF==mod_gene) %>% slice_max(n=1, importance) %>%
#  mutate(target=mod_gene)
#imp.df <- bind_rows(o2, o1) %>%
#  mutate(y_lab=factor(target, levels=c(rev(o1$target), mod_gene)),
#         bar_format= factor(ifelse(target==mod_gene, "self", "normal"), levels=c("self","normal")))
#
#plot.genes = imp.df$target
#length(plot.genes) 

p3 <- getDotplot(plot.genes, de.df)
p3.1 <- getMeanRatioBar(plot.genes, sce_summ)
#p3.2 <- getDetectedBoxplot(plot.genes, spe_summ)
#p3.3 <- ggplot(mutate(imp.df, y_lab=factor(y_lab, levels=rev(plot.genes))), 
#               aes(y=y_lab, x=importance, fill=bar_format, lty=bar_format, color=bar_format))+
#  geom_bar(stat="identity", linewidth=.3)+
#  scale_fill_manual(values=c("white","grey"), guide="none")+
#  scale_linetype_manual(values=c(2,1), guide="none")+
#  scale_color_manual(values=c("black","grey"), guide="none")+
#  labs(title=" ")+scale_y_discrete(position="right")+
#  theme_minimal()+theme(axis.text.y=element_blank(), 
#                        axis.title.y=element_blank(),
#                        axis.text.x=element_text(size=8), 
#                        panel.grid.minor=element_blank(), #panel.grid.major.y=element_blank(),
#                        plot.margin = margin(.2,.0,1.8,0, unit="cm"))

# mean ratio bar for InhN types
load("processed-data/publication/sce_SZBDMulti-seq_control_InhN_dotplot.Rdata")

ordered_genes = plot.genes
gids = rownames(sce_summ)[rowData(sce_summ)$gene_name %in% ordered_genes]
names(gids) = rowData(sce_summ)[gids,"gene_name"]	

df2 = tibble::rownames_to_column(as.data.frame(assay(sce_summ, "logcounts.mean")[gids,]))
colnames(df2) = c("gene_id", as.character(colData(sce_summ)$inhn_type))

bar.df = tidyr::pivot_longer(df2, all_of(levels(colData(sce_summ)$inhn_type)), 
                             names_to="cellType", values_to="mean.expr") %>%
	  left_join(as.data.frame(rowData(sce_summ)[gids,c("gene_id","gene_name")]))

bar.df$cellType = factor(bar.df$cellType, levels=levels(sce_summ$inhn_type))
bar.df$gene_name = factor(bar.df$gene_name, levels=rev(ordered_genes))

inhn.col.pal = c("CGE CNR1"="#5E646E","CGE LAMP5"="grey70","MGE PV"="#897d74","MGE SST"="#d6cac0")

p3.2 = ggplot(bar.df, aes(y=gene_name, x=mean.expr, fill=cellType))+
  geom_bar(stat="identity", position="fill")+
  scale_fill_manual(values=inhn.col.pal)+scale_y_discrete(position="right")+
  labs(title=" ", x="mean\nexpr")+guides(fill = guide_legend(nrow = 2, byrow = TRUE, position="bottom"))+
  theme_minimal()+theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
                        axis.title.y=element_blank(), panel.grid.minor=element_blank(),
                        #legend.position="bottom", 
			legend.title=element_blank(),
                        legend.key.size = unit(10,"pt"))

ggsave(file=paste0("plots/publication/Figure_nrn/",mod_gene,"_dotplot.pdf"), 
       arrangeGrob(grobs=list(p3, p3.1, p3.2), layout_matrix=matrix(c(1,1,1,1,1,1,2,3), ncol=8)),
       height=3.5, width=6.5)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
