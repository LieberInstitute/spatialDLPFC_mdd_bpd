formatManhattan <- function(.target_group, .target_sex, .orig_group, source,
                            gmt.list, gmtPadj, gsea_terms) {
  if(source=="smoothed") {
    match_group = names(sm_match)[sm_match==.orig_group]
    source2 = "sm"
    lf.df = filter(restr.results_sm, smoothed==.orig_group, group==.target_group, sex==.target_sex) %>% 
      group_by(gene_name) %>% summarise(avg_lf=mean(logFC)) %>% arrange(desc(avg_lf))
    lf.order = lf.df$avg_lf
    names(lf.order) = lf.df$gene_name
  }
  if(source=="seurat") {
    match_group = names(se_match)[se_match==.orig_group]
    source2 = "se"
    lf.df = filter(restr.results_se, seurat_label_f==.orig_group, group==.target_group, sex==.target_sex) %>%
      group_by(gene_name) %>% summarise(avg_lf=mean(logFC)) %>% arrange(desc(avg_lf))
    lf.order = lf.df$avg_lf
    names(lf.order) = lf.df$gene_name
  }
  sig.both.dir = c(pairList[[paste(.target_group, .target_sex, "dn", sep="_")]][[1]][[source]][[.orig_group]],
                   pairList[[paste(.target_group, .target_sex, "up", sep="_")]][[1]][[source]][[.orig_group]])
  paired.both.dir = c(pairList[[paste(.target_group, .target_sex, "dn", sep="_")]][[2]][[match_group]],
                      pairList[[paste(.target_group, .target_sex, "up", sep="_")]][[2]][[match_group]])
  
  results = gmtPadj[[paste(source2, .orig_group, sep="_")]]
  results = ungroup(results[results$pathway %in% gsea_terms,]) 
  
  plot.df = do.call(rbind, lapply(1:nrow(results), function(x) {
    term1 = results$pathway[x]
    as.data.frame(list("id"=names(lf.order),
                       "weight"=lf.order,
                       "rank"=1:length(lf.order),
                       "term"=term1,
                       "present"=names(lf.order) %in% gmt.list[[results[["pathway"]][x]]],
                       y=x-1, yend=x))
  })) %>%
    mutate(present2= ifelse(id %in% sig.both.dir & present==T, "LR", "NS"),
           present2= ifelse(id %in% paired.both.dir & present==T, "paired", present2),
           term = factor(term, levels=gsea_terms))
  
  
  label.df= as.data.frame(list(y_label= sapply(1:nrow(results), function(x) paste0("NES= ",round(results[["NES"]][x],2),
                                                                               " (",results[["size"]][x],"/ ",length(gmt.list[[results[["pathway"]][x]]]),")\n",
                                                                               "padj= ",format(results[["padj"]][x], scientific=T, digits=2))),
                               term=factor(results$pathway, levels=gsea_terms)))
  
  return(list("plotting"=plot.df,
              "labels"=label.df))

}
