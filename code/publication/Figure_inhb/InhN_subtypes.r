setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(edgeR)
	library(scater)
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
})
set.seed(123)

# norm psedobulk
load("processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo-InhN_indivID-inhn-type.Rdata")

rowData(sce_pseudo)$high_expr_group_individualID <- filterByExpr(sce_pseudo, group = sce_pseudo$individualID)
rowData(sce_pseudo)$high_expr_group_cluster <- filterByExpr(sce_pseudo, group = sce_pseudo$inhn_type)

sce_keep.genes = rownames(sce_pseudo)[rowData(sce_pseudo)$high_expr_group_individualID & rowData(sce_pseudo)$high_expr_group_cluster]

tmp = calcNormFactors(sce_pseudo)
x = cpm(tmp, log=T, prior.count=4)
dimnames(x) <- dimnames(sce_pseudo)
logcounts(sce_pseudo) <- x

sub_markers = c("ADARB2","CNR1","LAMP5","LHX6","SST","PVALB")#,"DLX1","ARX")
for(i in sub_markers) colData(sce_pseudo)[[i]] = logcounts(sce_pseudo)[rowData(sce_pseudo)$gene_name==i,]

cdata = as.data.frame(colData(sce_pseudo))
cdata2 = tidyr::pivot_longer(cdata, all_of(sub_markers), names_to="gene_name", values_to="logcounts")
cdata2$gene_name = factor(cdata2$gene_name, levels=sub_markers)

inhn.col.pal = c("CGE CNR1"="#5E646E","CGE LAMP5"="grey70","MGE PV"="#897d74","MGE SST"="#d6cac0")

p1 <- ggplot(cdata2, aes(x=inhn_type, y=logcounts, color=inhn_type))+
  ggbeeswarm::geom_quasirandom(size=.5)+
  facet_wrap(vars(gene_name))+ylim(0,15)+
  scale_color_manual(values=inhn.col.pal)+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), legend.position="none", text=element_text(size=6),
	panel.grid.minor=element_blank())

# azimuth collapse
load("processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_InhN-only.Rdata")

cdata = as.data.frame(colData(sce_con))
tmp = group_by(cdata, azimuth, inhn_type) %>% tally()

tmp$azimuth = factor(tmp$azimuth, levels=c("Vip","Sncg","Pax6","Lamp5","Lamp5 Lhx6","Sst","Sst Chodl","Pvalb","Chandelier"))

p2 <- ggplot(tmp, aes(x=azimuth, y=n, fill=inhn_type))+
  geom_bar(stat="identity", color="black")+labs(y="# nuclei")+
  scale_fill_manual(values=inhn.col.pal)+
  theme_bw()+theme(legend.position="none", text=element_text(size=6), panel.grid.minor=element_blank(),
	panel.grid.major.x=element_blank())

ggsave(file="plots/publication/Figure_nrn/InhN-subtypes_lineage-groups.pdf", 
	grid.arrange(p2, p1, layout_matrix=matrix(c(1,2,2))), width=3.5, height=5)

stop("First plot only")

# aggregate for mean ratio bar
sce_summ = aggregateAcrossCells(sce_con, ids=colData(sce_con)[,c("inhn_type")],
                            statistics=c("mean","prop.detected"),
                            use.assay.type="logcounts")

save(sce_summ, file="processed-data/publication/sce_SZBDMulti-seq_control_InhN_dotplot.Rdata")

# extra plot justifying ordered genes in violin form
ordered_genes = c("PVALB","RGS5","KCNS3","TAC1","TRBC2",
"GAD1","GAD2","ZNF385D","PNOC","LGI2","SLC32A1","SLC6A1",
"DLX6-AS1","RELN","VIP")

for(i in ordered_genes) colData(sce_pseudo)[[i]] = logcounts(sce_pseudo)[rowData(sce_pseudo)$gene_name==i,]

cdata = as.data.frame(colData(sce_pseudo))
cdata2 = tidyr::pivot_longer(cdata, all_of(gsub("DLX6-AS1","DLX6\\.AS1",ordered_genes)), names_to="gene_name", values_to="logcounts")
cdata2$gene_name = factor(cdata2$gene_name, levels=gsub("DLX6-AS1","DLX6\\.AS1",ordered_genes), labels=ordered_genes)

p3 <- ggplot(cdata2, aes(x=inhn_type, y=logcounts, color=inhn_type))+
  ggbeeswarm::geom_quasirandom()+
  facet_wrap(vars(gene_name))+ylim(0,15)+
  scale_color_manual(values=RColorBrewer::brewer.pal("Paired", n=4))+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), legend.position="bottom", text=element_text(size=6))

ggsave(file="plots/publication/Figure_nrn/InhN-subtypes_GAD1-module.pdf",
        p3, width=5, height=8)

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
