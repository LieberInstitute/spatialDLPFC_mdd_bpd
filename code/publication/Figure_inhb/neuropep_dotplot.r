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

plot.genes = c("CRH","SST","CORT")

p3 <- getDotplot(plot.genes, de.df2)+theme(axis.text.y=element_text(size=7))
p3.1 <- getMeanRatioBar(plot.genes, sce_summ)


## mean ratio bar for InhN types
#load("processed-data/publication/sce_SZBDMulti-seq_control_InhN_dotplot.Rdata")
#
#ordered_genes = plot.genes
#gids = rownames(sce_summ)[rowData(sce_summ)$gene_name %in% ordered_genes]
#names(gids) = rowData(sce_summ)[gids,"gene_name"]	
#
#df2 = tibble::rownames_to_column(as.data.frame(assay(sce_summ, "logcounts.mean")[gids,]))
#colnames(df2) = c("gene_id", as.character(colData(sce_summ)$inhn_type))
#
#bar.df = tidyr::pivot_longer(df2, all_of(levels(colData(sce_summ)$inhn_type)), 
#                             names_to="cellType", values_to="mean.expr") %>%
#	  left_join(as.data.frame(rowData(sce_summ)[gids,c("gene_id","gene_name")]))
#
#bar.df$cellType = factor(bar.df$cellType, levels=levels(sce_summ$inhn_type))
#bar.df$gene_name = factor(bar.df$gene_name, levels=rev(ordered_genes))
#
#inhn.col.pal = c("CGE CNR1"="#5E646E","CGE LAMP5"="grey70","MGE PV"="#897d74","MGE SST"="#d6cac0")
#
#p3.2 = ggplot(bar.df, aes(y=gene_name, x=mean.expr, fill=cellType))+
#  geom_bar(stat="identity", position="fill")+
#  scale_fill_manual(values=inhn.col.pal)+scale_y_discrete(position="right")+
#  labs(title=" ", x="mean\nexpr")+guides(fill = guide_legend(nrow = 2, byrow = TRUE, position="bottom"))+
#  theme_minimal()+theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
#                        axis.title.y=element_blank(), panel.grid.minor=element_blank(),
#                        #legend.position="bottom", 
#			legend.title=element_blank(),
#                        legend.key.size = unit(10,"pt"))

ggsave(file="plots/publication/Figure_nrn/neuropep_dotplot.pdf", 
       arrangeGrob(grobs=list(p3, p3.1), layout_matrix=matrix(c(1,1,1,1,1,2), ncol=6)),
       height=2, width=3.5)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
