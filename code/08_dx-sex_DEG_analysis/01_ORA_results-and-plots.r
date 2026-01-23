setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(clusterProfiler)
	library(enrichR)
	library(pheatmap)
	library(igraph)
	library(gridExtra)
	library(grid)
	library(gtable)
})
set.seed(123)

source("code/08_dx-sex_DEG_analysis/custom_functions.r")

results_set = "smoothed-k9-1663"
#results_set = "seurat-pc30"

# format gmt for enricher
setEnrichrSite("Enrichr") # Human genes
dbs <- listEnrichrDbs()

gmt = .read_gmt("Reactome_2022")
gmt.id.split = " R-HSA-"

## separate term and RSA for TERM2NAME
term2name = do.call(rbind.data.frame, strsplit(as.character(gmt$term), split=gmt.id.split))
colnames(term2name) <- c("Term","ID")
term2name$ID = paste0(substr(gmt.id.split, start=2, stop=nchar(gmt.id.split)), term2name$ID)
term.id.long = term2name$ID
term2name = distinct(term2name[,c("ID","Term")])

## dframe of goID and gene name for TERM2GENE
term2gene = cbind.data.frame("ID"=term.id.long, "geneID"=gmt$gene)


# L-A ORA
la.degs = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))
#la.degs_t = filter(la.degs, n_ttest_sig>0)

gene_universe = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,
                                "_rev-gene-input_F-test.csv"))$gene_name

ora = enricher(la.degs$gene_name, universe=gene_universe, 
               TERM2GENE = term2gene,
               TERM2NAME = term2name)
ora@result <- ora@result[ora@result$p.adjust<.05,]

saveRDS(ora, paste0("processed-data/08_dx-sex_DEG_analysis/layer-adjusted_", results_set,"_F-test-padj05_Reactome-ORA.rda"))
cat("\nReactome ORA results for L-A saved to:", paste0("processed-data/08_dx-sex_DEG_analysis/layer-adjusted_", 
  results_set,"_F-test-padj05_Reactome-ORA.rda"), "\n\n")

# L-R ORA
lr.degs = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))
#lr.degs_t = filter(lr.degs, n_ttest_sig>0)

ora_lr = enricher(lr.degs$gene_name, universe=gene_universe, 
               TERM2GENE = term2gene,
               TERM2NAME = term2name)
ora_lr@result <- ora_lr@result[ora_lr@result$p.adjust<.05,]

saveRDS(ora_lr, paste0("processed-data/08_dx-sex_DEG_analysis/layer-restricted_", results_set,"_F-test-padj05_Reactome-ORA.rda"))
cat("\nReactome ORA results for L-R saved to:", paste0("processed-data/08_dx-sex_DEG_analysis/layer-restricted_", 
  results_set,"_F-test-padj05_Reactome-ORA.rda"), "\n\n")


# L-A heatmap
test_set = ora@result
p.names = sapply(test_set$Description, function(x) {
  c1 = unlist(strwrap(x, width=75))
  if(length(c1)>1) {
    return(paste(c1[[1]],"[...]"))
  } else {
    return(c1[[1]])
  }
})
test_set$name <- p.names
expand_genes = tidyr::separate_rows(test_set, geneID, sep="/")

la.degs_t = filter(la.degs, n_ttest_sig>0)
tsig = la.degs_t$gene_name

tsig_genes = filter(expand_genes, geneID %in% tsig)
df1 = tidyr::pivot_wider(tsig_genes[,c("name", "FoldEnrichment", "geneID")] , 
                          names_from="geneID", values_from="FoldEnrichment", values_fill=0)
m1 = as.matrix(df1[,-1])
rownames(m1) <- df1$name

#color limits for NES
color.limits = ceiling(max(abs(m1)))
col1 = colorRampPalette(RColorBrewer::brewer.pal(n=5, "OrRd"))(color.limits)
col1[1] = "grey90"

phm1 = pheatmap(m1, color=col1,
                breaks= seq(0, color.limits),
                treeheight_col = 12, treeheight_row = 12, border_color = "grey",
                silent=T, 
                fontsize_row=7, fontsize_col=7, angle_col=90,
                main="Fill: FoldEnrichment (0 if not sig.)")


# L-R heatmap
test_set = ora_lr@result
p.names = sapply(test_set$Description, function(x) {
  c1 = unlist(strwrap(x, width=75))
  if(length(c1)>1) {
    return(paste(c1[[1]],"[...]"))
  } else {
    return(c1[[1]])
  }
})
test_set$name <- p.names
expand_genes = tidyr::separate_rows(test_set, geneID, sep="/")

lr.degs_t = filter(lr.degs, n_ttest_sig>0)
tsig = lr.degs_t$gene_name

tsig_genes2 = filter(expand_genes, geneID %in% tsig)
df2 = tidyr::pivot_wider(tsig_genes2[,c("name", "FoldEnrichment", "geneID")] , 
                          names_from="geneID", values_from="FoldEnrichment", values_fill=0)
m2 = as.matrix(df2[,-1])
rownames(m2) <- df2$name

## color limits for NES
color.limits = ceiling(max(abs(m2)))
col1 = colorRampPalette(RColorBrewer::brewer.pal(n=5, "OrRd"))(color.limits)
col1[1] = "grey90"

phm2 = pheatmap(m2, color=col1,
                breaks= seq(0, color.limits),
                treeheight_col = 12, treeheight_row = 12, border_color = "grey",
                silent=T, 
                fontsize_row=7, fontsize_col=7, angle_col=90,
                main="Fill: FoldEnrichment (0 if not sig.)")

# make term2gene igraphs for all dx*sex groups for LA
comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")
names(comparisons) <- comparisons

igraphList <- lapply(comparisons, function(x) {
  coef1 = la.degs_t[,c("gene_id","gene_name","med_prop.spots.detected",
                       grep(x, colnames(la.degs_t), value=T))]
  colnames(coef1)[4:5] = c("logFC","ttest")
  decreased = filter(coef1, logFC<0, ttest!="NS")$gene_name
  increased = filter(coef1, logFC>0, ttest!="NS")$gene_name
  
  #for grid table with number of genes
  ngenes = c("n_dn"=length(decreased), "n_up"=length(increased),
             "n_dn_annot"=length(intersect(decreased, tsig_genes$geneID)),
             "n_up_annot"=length(intersect(increased, tsig_genes$geneID)))
  
  tmp = bind_rows(filter(tsig_genes, geneID %in% increased) %>% mutate(NES=FoldEnrichment),
                  filter(tsig_genes, geneID %in% decreased) %>% mutate(NES=-FoldEnrichment)) %>%
    select(term=ID, name, NES, gene_name=geneID)

  tmp = left_join(tmp, la.degs_t[,c("gene_name","med_prop.spots.detected")])
  tmp$node_size = ifelse(tmp$med_prop.spots.detected>.05, 6, 3)
  
    outlist <- term2GeneIGRAPH(tmp, source=FALSE, node_size="node_size", layout_style="kk", 
                             text_title=paste0("L-A DEGs (",results_set,"): ", x), color_by="gene") 
  
  return(list("ngenes"=ngenes, "igraphL"=outlist))
})


# summarise L-A DEG results in gtables
mytheme <- ttheme_default(
  core = list(fg_params=list(cex = 1)),
  colhead = list(fg_params=list(cex = 1)))

grobList <- lapply(comparisons, function(x) {
  tmp = unlist(strsplit(x, "_"))
  target_group = tmp[[1]]
  target_sex = tmp[[2]]
  
  t1 = data.frame("set"=c("# L-A DEGs","# L-A DEGs annot.\nby F-test ORA"),
                  "dn"=igraphList[[x]]$ngenes[c("n_dn","n_dn_annot")],
                  "up"=igraphList[[x]]$ngenes[c("n_up","n_up_annot")]
                  )
  gt1 = tableGrob(t1, rows = NULL, cols=c("","Decreased","Increased"),
                  theme = mytheme)
  
  #add title
  title <- textGrob(paste(target_group, target_sex),gp=gpar(fontsize=20))
  padding <- unit(5,"mm")
  table <- gtable_add_rows(
    gt1, 
    heights = grobHeight(title) + padding,
    pos = 0)
  table <- gtable_add_grob(
    table, 
    title, 
    1, 1, 1, ncol(table))
  return(table)
})


# save plots to PDF
pdf(file=paste0("plots/08_dx-sex_DEG_analysis/Reactome-ORA_", results_set, 
                "_F-test-padj05.pdf"), 
    width=12, height=12)
grid.arrange(phm1[[4]], top=paste0("Layer-adjusted (", results_set, ")"))
grid.arrange(grobList[[1]], grobList[[2]],
             grobList[[3]], grobList[[4]],
             grobList[[5]], grobList[[6]], 
             layout_matrix=matrix(1:6, ncol=2, byrow=T),
	top=paste0("Layer-adjusted (", results_set, ")"),
	bottom="DEGs: F test adj. p<.05 and moderated t-test adj.p <.05"
)
for(i in igraphList) {
  plot(i$igraphL$igraph, layout=i$igraphL$layout, vertex.label.family="sans", #vertex.size=6, 
	edges.curved=T,
       vertex.frame.width=0,
       main=i$igraphL[["title_text"]],
       sub=paste0("Reactome annotations determined by ORA with n= ", nrow(la.degs)," genes with F-test adj. p<.05",
	"\nSmaller genes have median zero proportion >=5% of spots."))
  legend("bottomleft", legend=c("Decreased","Increased","Term"),
         pch=16, pt.cex=1, cex=.7,
         col=c("#CFEBF7","#FFC0B5","grey85"))
}
grid.arrange(phm2[[4]], top=paste0("Layer-restricted (", results_set, ")"))
dev.off()

cat("\n\nReactome ORA results plots saved to:", paste0("plots/08_dx-sex_DEG_analysis/Reactome-ORA_", 
  results_set,  "_F-test-padj05.pdf"),"\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
