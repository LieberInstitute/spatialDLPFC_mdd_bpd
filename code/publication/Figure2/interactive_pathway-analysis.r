library(dplyr)
library(clusterProfiler)
#library(enrichR)
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


test1 = group_by(expand_genes3, geneID, model) %>% tally(name="n_terms") %>%
  tidyr::pivot_wider(names_from="model", values_from="n_terms", values_fill=0)
test1$is_LA_F.sig = test1$geneID %in% la.degs$gene_name
test1$is_LR_F.sig = test1$geneID %in% lr.degs$gene_name
test1$F_set = factor(paste(test1$is_LA_F.sig, test1$is_LR_F.sig), levels=c("TRUE FALSE","FALSE TRUE","TRUE TRUE"),
                     labels=c("L-A only","L-R only","both"))
table(test1$F_set)
#L-A only L-R only     both 
#132       26       65


test2 = distinct(expand_genes3, ID, Description, name, model, Count) %>%
  tidyr::pivot_wider(names_from="model", values_from="Count") %>%
  mutate(is_LA_term= !is.na(LA), is_LR_term= !is.na(LR), 
         F_set=factor(paste(is_LA_term, is_LR_term), levels=c("TRUE FALSE","FALSE TRUE","TRUE TRUE"),
                      labels=c("L-A only","L-R only","both")))
table(test2$F_set)
#L-A only L-R only     both 
#25       33       30 


test3 = group_by(expand_genes3, ID, Description, name, geneID) %>% 
  summarise(avg_FE=mean(FoldEnrichment), avg_Count=mean(Count),
            avg_zScore=mean(zScore))
ggplot(distinct(test3, ID, Description, avg_FE, avg_Count), aes(x=avg_Count, y=avg_FE))+
  geom_point()




df = tidyr::pivot_wider(test3[,c("name", "avg_FE", "geneID")], 
                         names_from="geneID", values_from="avg_FE", values_fill=0)
m1 = as.matrix(df[,-1])
rownames(m1) <- df$name

## color limits for NES
color.limits = ceiling(max(abs(m1)))
col1 = colorRampPalette(RColorBrewer::brewer.pal(n=5, "OrRd"))(color.limits)
col1[1] = "white"

annot_row = data.frame("term_in"=test2$F_set, row.names=test2$name)
annot_col = data.frame("gene_in"=test1$F_set, row.names=test1$geneID)

annot_colors= list("term_in"=c("L-A only"="gold", "L-R only"="pink2", "both"="grey"),
                   "gene_in"=c("L-A only"="gold", "L-R only"="pink2", "both"="grey"))
phm1 = pheatmap(m1, color=col1,
                breaks= seq(0, color.limits), 
                method="average",clustering_distance_rows="binary", clustering_distance_cols="binary",
                treeheight_col = 12, treeheight_row = 12, 
                #silent=T, 
                show_colnames=F,
                fontsize_row=7,
                annotation_row=annot_row, annotation_col=annot_col, annotation_colors=annot_colors,
                annotation_names_row=F, annotation_names_col=F,
                main="Fill: FoldEnrichment (0 if not sig.)")




# add L-A GO:BP
ora.bp = readRDS(paste0("processed-data/08_dx-sex_DEG_analysis/layer-adjusted_", results_set,"_F-test-padj05_GO-BP-ORA.rda"))
#ora_lr = readRDS(paste0("processed-data/08_dx-sex_DEG_analysis/layer-restricted_", results_set,"_F-test-padj05_GO-BP-ORA.rda"))
## L-R only has a few sig terms so not included

# L-A heatmap
bp_set = ora.bp@result
p.names = sapply(bp_set$Description, function(x) {
  c1 = unlist(strwrap(x, width=75))
  if(length(c1)>1) {
    return(paste(c1[[1]],"[...]"))
  } else {
    return(c1[[1]])
  }
})
bp_set$name <- p.names
expand_genes4 = tidyr::separate_rows(bp_set, geneID, sep="/")



expand_genes5 = bind_rows(mutate(expand_genes4, model="LA_BP"),
                          expand_genes3)


test4 = group_by(expand_genes5, geneID, model) %>% tally(name="n_terms") %>%
  tidyr::pivot_wider(names_from="model", values_from="n_terms", values_fill=0)
test4$is_LA_F.sig = test4$geneID %in% la.degs$gene_name
test4$is_LR_F.sig = test4$geneID %in% lr.degs$gene_name
test4$F_set = factor(paste(test4$is_LA_F.sig, test4$is_LR_F.sig), levels=c("TRUE FALSE","FALSE TRUE","TRUE TRUE"),
                     labels=c("L-A only","L-R only","both"))
table(test4$F_set)
#L-A only L-R only     both 
#209       26       98 


test5 = distinct(expand_genes5, ID, Description, name, model, Count) %>%
  tidyr::pivot_wider(names_from="model", values_from="Count") %>%
  mutate(is_LA_term= !is.na(LA), is_LR_term= !is.na(LR), is_BP_term= !is.na(LA_BP),
         F_set=factor(paste(is_BP_term, is_LA_term, is_LR_term), 
                      levels=c("TRUE FALSE FALSE",
                               "FALSE TRUE FALSE","FALSE FALSE TRUE","FALSE TRUE TRUE"),
                      labels=c("BP L-A",
                               "Rctme L-A only","Rctme L-R only","Rctme both")))
table(test5$F_set)
#BP L-A Rctme L-A only Rctme L-R only     Rctme both 
#.   77             25             33             30

group_by(test5, name) %>% tally() %>% filter(n>1)
#"Translation is there twice so just change the GO:BP version
test5[test5$name=="Translation" & test5$is_BP_term, "name"] = "Translation (GO-BP)"

test6 = group_by(expand_genes5, ID, Description, name, geneID) %>% 
  summarise(avg_FE=mean(FoldEnrichment), avg_Count=mean(Count),
            avg_zScore=mean(zScore))
change.names = intersect(grep("GO", test6$ID), grep("^Translation$", test6$name))
test6[change.names,"name"] = "Translation (GO-BP)"
ggplot(distinct(test6, ID, Description, avg_FE, avg_Count), aes(x=avg_Count, y=avg_FE))+
  geom_point()

df2 = tidyr::pivot_wider(test6[,c("name", "avg_FE", "geneID")], 
                        names_from="geneID", values_from="avg_FE", values_fill=0)
m2 = as.matrix(df2[,-1])
rownames(m2) <- df2$name

## color limits for NES
color.limits = ceiling(max(abs(m2)))
col1 = colorRampPalette(RColorBrewer::brewer.pal(n=5, "OrRd"))(color.limits)
col1[1] = "white"

annot_row = data.frame("term_in"=test5$F_set, row.names=test5$name)
annot_col = data.frame("gene_in"=test4$F_set, row.names=test4$geneID)

annot_colors= list("term_in"=c("BP L-A"="palegoldenrod","Rctme L-A only"="gold", "Rctme L-R only"="pink2", "Rctme both"="grey"),
                   "gene_in"=c("L-A only"="gold", "L-R only"="pink2", "both"="grey"))
phm2 = pheatmap(m2, color=col1,
                breaks= seq(0, color.limits), 
                method="average",clustering_distance_rows="binary", clustering_distance_cols="binary",
                treeheight_col = 12, treeheight_row = 12, 
                #silent=T, 
                show_colnames=F,
                fontsize_row=7,
                annotation_row=annot_row, annotation_col=annot_col, annotation_colors=annot_colors,
                annotation_names_row=F, annotation_names_col=F,
                main="Fill: FoldEnrichment (0 if not sig.)")

phm3 = pheatmap(m2, color=col1,
                breaks= seq(0, color.limits), 
                method="average",clustering_distance_rows="binary", clustering_distance_cols="binary",
                treeheight_col = 12, treeheight_row = 12, 
                silent=T, 
                show_colnames=F, show_rownames=F,
                #fontsize_row=7,
                annotation_row=annot_row, annotation_col=annot_col, annotation_colors=annot_colors,
                annotation_names_row=F, annotation_names_col=F,
                main="Fill: FoldEnrichment (0 if not sig.)")

test7 = bind_rows(filter(expand_genes5, model %in% c("LA", "LA_BP"), geneID %in% filter(la.degs, n_ttest_sig>0)$gene_name),
                  filter(expand_genes5, model=="LR", geneID %in% filter(lr.degs, n_ttest_sig>0)$gene_name)) %>%
  group_by(ID, Description, name, geneID) %>% 
  summarise(avg_FE=mean(FoldEnrichment), avg_Count=mean(Count))
change.names = intersect(grep("GO", test7$ID), grep("^Translation$", test7$name))
test7[change.names,"name"] = "Translation (GO-BP)"
#ggplot(distinct(test6, ID, Description, avg_FE, avg_Count), aes(x=avg_Count, y=avg_FE))+
#  geom_point()

df3 = tidyr::pivot_wider(test7[,c("name", "avg_FE", "geneID")], 
                         names_from="geneID", values_from="avg_FE", values_fill=0)
m3 = as.matrix(df3[,-1])
rownames(m3) <- df3$name

phm4 = pheatmap(m3, color=col1,
                breaks= seq(0, color.limits), 
                method="average",clustering_distance_rows="binary", clustering_distance_cols="binary",
                treeheight_col = 12, treeheight_row = 12, 
                silent=T, 
                show_colnames=F,
                fontsize_row=7,
                annotation_row=annot_row, annotation_col=annot_col, annotation_colors=annot_colors,
                annotation_names_row=F, annotation_names_col=F,
                main="Fill: FoldEnrichment (0 if not sig.)")

phm5 = pheatmap(m3, color=col1,
                breaks= seq(0, color.limits), 
                method="average",clustering_distance_rows="binary", clustering_distance_cols="binary",
                treeheight_col = 12, treeheight_row = 12, 
                silent=T, 
                show_colnames=F, show_rownames=F,
                annotation_row=annot_row, annotation_col=annot_col, annotation_colors=annot_colors,
                annotation_names_row=F, annotation_names_col=F,
                main="Fill: FoldEnrichment (0 if not sig.)")


pdf(file="plots/publication/Figure2/pathway-analysis_ORA_heatmap.pdf", height=12, width=8)
plot(phm2[[4]])
plot(phm3[[4]])
plot(phm4[[4]])
plot(phm5[[4]])
dev.off()

