setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(dplyr)
  library(clusterProfiler)
  library(pheatmap)
})
set.seed(123)

results_set = "smoothed-k9-1663"
la.degs = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))
lr.degs = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))

ora = readRDS(paste0("processed-data/08_dx-sex_DEG_analysis/layer-adjusted_", results_set,"_F-test-padj05_Reactome-ORA.rda"))
ora_lr = readRDS(paste0("processed-data/08_dx-sex_DEG_analysis/layer-restricted_", results_set,"_F-test-padj05_Reactome-ORA.rda"))

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


# L-R heatmap
test_set2 = ora_lr@result
p.names = sapply(test_set2$Description, function(x) {
  c1 = unlist(strwrap(x, width=75))
  if(length(c1)>1) {
    return(paste(c1[[1]],"[...]"))
  } else {
    return(c1[[1]])
  }
})
test_set2$name <- p.names
expand_genes2 = tidyr::separate_rows(test_set2, geneID, sep="/")

expand_genes3 = bind_rows(mutate(expand_genes, model="LA"),
                          mutate(expand_genes2, model="LR"))


tmp = filter(expand_genes3, model=="LA")
tmp2 = tidyr::pivot_wider(tmp[,c("ID","name","FoldEnrichment","geneID")], names_from="geneID", values_from="FoldEnrichment", values_fill=0)
mtx1 = as.matrix(tmp2[,3:ncol(tmp2)])
rownames(mtx1) <- tmp2$name
phm1 <- pheatmap(mtx1, fontsize=6, main="Layer-adjusted F-test ORA: Reactome")

tt = filter(la.degs, gene_name %in% colnames(mtx1))

#col_annot = as.data.frame(tt[,c(18,20,22,19,21,23)], row.names = tt$gene_name)
#annot_colors= lapply(colnames(col_annot), function(x) c("NS"="white","padj<.05"="black","padj<.01"="black","padj<.0001"="black"))
#names(annot_colors) = colnames(col_annot)
#phm1 <- pheatmap(mtx1, fontsize=6, main="Layer-adjusted F-test ORA: Reactome",
#                 color = colorRampPalette(c("white","black"))(100),
#                 annotation_col=col_annot[,ncol(col_annot):1], annotation_colors=annot_colors, 
#                 annotation_legend = F, angle_col=90, treeheight_row = 20, treeheight_col = 20)


col_annot_alt = as.data.frame(tt[,c(1,3,5,2,4,6)], row.names = tt$gene_name)
summary(col_annot_alt)
#i will need to provide two color values for gradent for each coef and I will need to pick those colors to represent to min/max for each column based on a fixed reference that spans the whole spectrum/range

col_annot_alt = as.data.frame(signif(plyr::round_any(as.matrix(tt[,c(1,3,5,2,4,6)]), .05, f = round),2), row.names=tt$gene_name)
col_annot_alt = as.data.frame(apply(col_annot_alt, MARGIN=2, as.character), row.names=rownames(col_annot_alt))
# -3 to 3
minmax = colorRampPalette(rev(RColorBrewer::brewer.pal(n=7,"RdBu")))(120)
minmax_num = as.character(round(seq(-3, 3, length=121)[-1], 2))
names(minmax) = minmax_num

annot_colors_alt = lapply(colnames(col_annot_alt), function(x) minmax[col_annot_alt[,x]])
names(annot_colors_alt) = colnames(col_annot_alt)


phm1 = pheatmap(mtx1, fontsize=6, main="Layer-adjusted F-test ORA: Reactome",
         color = colorRampPalette(c("white","black"))(100), breaks=seq(0,12.2,length=101),
         annotation_col=col_annot_alt[,ncol(col_annot_alt):1], annotation_colors=annot_colors_alt, 
         annotation_legend = F, angle_col=90, treeheight_row = 20, treeheight_col = 20)

phm1_norow = pheatmap(mtx1, fontsize=6, main="Layer-adjusted F-test ORA: Reactome",
                color = colorRampPalette(c("white","black"))(100), breaks=seq(0,12.2,length=101),
                annotation_col=col_annot_alt[,ncol(col_annot_alt):1], annotation_colors=annot_colors_alt, 
                annotation_legend = F, angle_col=90, treeheight_row = 20, treeheight_col = 20,
                show_rownames = F, fontsize_col=4)

fake_mtx = as.matrix(tt[,c(1,3,5,2,4,6)], row.names=tt$gene_name)
phm1_legend = pheatmap(fake_mtx, colors=colorRampPalette(rev(RColorBrewer::brewer.pal(n=7,"RdBu")))(100),
         breaks=seq(-3,3,length=101), main="Fake mtx for legend (logFC values)")

tmplr = filter(expand_genes3, model=="LR")
tmplr2 = tidyr::pivot_wider(tmplr[,c("ID","name","FoldEnrichment","geneID")], names_from="geneID", values_from="FoldEnrichment", values_fill=0)
mtx2 = as.matrix(tmplr2[,3:ncol(tmplr2)])
rownames(mtx2) <- tmplr2$name
phm2 <- pheatmap(mtx2, fontsize=6, main="Layer-restricted F-test ORA: Reactome")

colnames(lr.degs)
tt2 = filter(lr.degs, gene_name %in% colnames(mtx2))
col_annot2 = as.data.frame(tt2[,53:58]>0, row.names=tt2$gene_name)
col_annot2 = as.data.frame(apply(col_annot2, MARGIN=2, as.character), row.names=tt2$gene_name)
colnames(col_annot2) = c("L1","L2","L3.4","L5","L6","WM")

cpList <- readRDS("plots/colorPalettes.rds")
annot_colors2 = lapply(colnames(col_annot2), function(x) c("FALSE"="white", "TRUE"=cpList$smoothed.bright[[x]]))
names(annot_colors2) = colnames(col_annot2)

phm2 <- pheatmap(mtx2, fontsize=6, main="Layer-restricted F-test ORA: Reactome",
                 color = colorRampPalette(c("white","black"))(100), breaks=seq(0,12.2,length=101),
                 annotation_col=col_annot2[,ncol(col_annot2):1], annotation_colors=annot_colors2,
                 annotation_legend=F, angle_col=90, treeheight_row = 20, treeheight_col = 20)

phm2_norow = pheatmap(mtx2, fontsize=6, main="Layer-restricted F-test ORA: Reactome",
                      color = colorRampPalette(c("white","black"))(100), breaks=seq(0,12.2,length=101),
                      annotation_col=col_annot2[,ncol(col_annot2):1], annotation_colors=annot_colors2,
                      annotation_legend=F, angle_col=90, treeheight_row = 20, treeheight_col = 20,
                      show_rownames = F)

pdf(file="plots/publication/Figure2/pathway-analysis_ORA_heatmaps.pdf", height=6, width=6.5)
plot(phm1[[4]])
plot(phm2[[4]])
plot(phm1_norow[[4]])
plot(phm2_norow[[4]])
plot(phm1_legend[[4]])
dev.off()


## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()

